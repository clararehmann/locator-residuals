"""
run vector correlation between empirical datasets N times,
run permutation test for correlation N times
"""

import pandas as pd, numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.neighbors import NearestNeighbors
from scipy import stats
from scipy.stats import chi2
from scipy.stats import ttest_ind
from numpy.random import default_rng
from scipy.spatial import KDTree
import math, os

import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--ag', help='path to anopheles centroids')
parser.add_argument('--pf', help='path to plasmodium centroids')
parser.add_argument('--N', type=int, default=1000)
parser.add_argument('--out', help='path to output (appended with _correlation.txt and _permutation.txt')
args=parser.parse_args()

# functions
def veccor(u1, v1, u2, v2):
    """
    input is two vector timeseries
        u1: dx1 values
        v1: dy1 values
        u2: dx2 values
        v2: dy2 values
    cartesian coords
    """
    # stuff vectors into matrix
    x = np.array([u1,v1,u2,v2])
    # compute covariance matrix
    sigma = np.cov(x)
    # f terms are those that appear from left to right in Crosby et al
    f1 = sigma[0,0] * ((sigma[2,2] * sigma[1,3]**2) + (sigma[3,3] * sigma[1,2]**2))
    f2 = sigma[1,1] * ((sigma[2,2] * sigma[0,3]**2) + (sigma[3,3] * sigma[0,2]**2))
    f3 = 2 * sigma[0,1] * sigma[0,3] * sigma[1,2] * sigma[2,3]
    f4 = 2 * sigma[0,1] * sigma[0,2] * sigma[1,3] * sigma[2,3]
    f5 = -2 * sigma[0,0] * sigma[1,2] * sigma[1,3] * sigma[2,3]
    f6 = -2 * sigma[1,1] * sigma[0,2] * sigma[0,3] * sigma[2,3]
    f7 = -2 * sigma[2,2] * sigma[0,1] * sigma[0,3] * sigma[1,3]
    f8 = -2 * sigma[3,3] * sigma[0,1] * sigma[0,2] * sigma[1,2]
    # terms in denominator
    g1 = sigma[0,0] * sigma[1,1] - sigma[0,1]**2
    g2 = sigma[2,2] * sigma[3,3] - sigma[2,3]**2

    r = (f1 + f2 + f3 + f4 +  f5 + f6 + f7 + f8) / (g1 * g2)
    return r

def distance_matrix(df1, df2, max_distance=10):
    """
    compute pairwise distance matrix between two spatial datasets
    (for pairing samples)
    """
    # stack xy values
    xy1 = np.stack((df1.x, df1.y), axis=1)
    xy2 = np.stack((df2.x, df2.y), axis=1)
    # KDTrees
    tree1 = KDTree(xy1)
    tree2 = KDTree(xy2)
    # compute distance matrix
    distmat = tree1.sparse_distance_matrix(tree2, max_distance)
    distmat = distmat.toarray()
    return distmat

def find_pairs(df1, df2, distmat):
    """
    pair points across two datasets by nearest neighbors
    """
    # copy distance matrix so it doesn't get messed with
    dists = np.copy(distmat)

    # save indices
    idxs1 = []
    idxs2 = []

    # randomly shuffle who to pair
    for index1 in np.random.choice(np.arange(len(df1)), size=len(df1), replace=False):
        # distances from df2 points
        darr = dists[index1]
        # if any df2 neighbors,
        if sum(darr) > 0:
            # get the indices of closest
            index2 = np.where( darr == np.min(darr[np.nonzero(darr)]))[0]
            # randomly choose one
            index2 = np.random.choice(index2)
            # set it as 0 so it's not chosen again
            dists[:, index2] = 0

            # save to array
            idxs1.append(index1)
            idxs2.append(index2)
    return idxs1, idxs2

def get_vectors(df1, df2, idxs1, idxs2):
    """
    get neighboring points' dx dy vectors for veccor()
    """
    u1 = [] # dx1
    v1 = [] # dy1
    u2 = [] # dx2
    v2 = [] # dy2

    for i in range(len(idxs1)):
        idx1 = idxs1[i]
        idx2 = idxs2[i]
        u1.append(df1.loc[idx1, 'gc_x'] - df1.loc[idx1, 'x'])
        v1.append(df1.loc[idx1, 'gc_y'] - df1.loc[idx1, 'y'])
        u2.append(df2.loc[idx2, 'gc_x'] - df2.loc[idx2, 'x'])
        v2.append(df2.loc[idx2, 'gc_y'] - df2.loc[idx2, 'y'])

    return u1, v1, u2, v2

def get_angles(df1, df2, idxs1, idxs2):
    """
    get neighboring points' vector angles
    """
    a1 = []
    a2 = []

    for i in range(len(idxs1)):
        idx1 = idxs1[i]
        idx2 = idxs2[i]

        a1.append(math.degrees(math.atan2((df1.loc[idx1, 'gc_x'] - df1.loc[idx1, 'x']), (df1.loc[idx1, 'gc_y'] - df1.loc[idx1, 'y']))))
        a2.append(math.degrees(math.atan2((df2.loc[idx2, 'gc_x'] - df2.loc[idx2, 'x']), (df2.loc[idx2, 'gc_y'] - df2.loc[idx2, 'y']))))

    return a1, a2

def get_mags(df1, df2, idxs1, idxs2):
    """
    get neighboring points' vector magnitudes
    """
    w1 = []
    w2 = []

    for i in range(len(idxs1)):
        idx1 = idxs1[i]
        idx2 = idxs2[i]

        w1.append(np.hypot((df1.loc[idx1, 'gc_x'] - df1.loc[idx1, 'x']), (df1.loc[idx1, 'gc_y'] - df1.loc[idx1, 'y'])))
        w2.append(np.hypot((df2.loc[idx2, 'gc_x'] - df2.loc[idx2, 'x']), (df2.loc[idx2, 'gc_y'] - df2.loc[idx2, 'y'])))

    return w1, w2

# read in data
ag = pd.read_csv(args.ag, sep='\t')
pf = pd.read_csv(args.pf, sep='\t')
# distance matrix
distmat = distance_matrix(ag, pf)
 
# run correlation tests, shuffling within-region pairs
vec_correlation = np.empty(args.N)
ang_correlation = np.empty(args.N)
mag_correlation = np.empty(args.N)
 
for P in range(args.N):
    # neighbors
    agi, pfi = find_pairs(ag, pf, distmat)

    # data for pairs
    u1, v1, u2, v2 = get_vectors(ag, pf, agi, pfi)
    a1, a2 = get_angles(ag, pf, agi, pfi)
    w1, w2 = get_mags(ag, pf, agi, pfi)

    # test correlation
    vc = veccor(u1, v1, u2, v2)
    ac = stats.spearmanr(a1,a2)[0]
    wc = stats.spearmanr(w1,w2)[0]

    vec_correlation[P] = vc
    ang_correlation[P] = ac
    mag_correlation[P] = wc

# shove into dataframe and save
correlation = pd.DataFrame({'vector':vec_correlation, 'angle':ang_correlation, 'mag':mag_correlation})
correlation.to_csv(args.out+'_correlation.txt', sep='\t')

# run permutation tests, shuffling all data pairs
vec_correlation = np.empty(args.N)
ang_correlation = np.empty(args.N)
mag_correlation = np.empty(args.N)

for P in range(args.N):
    # randomize pairs of points (ag < pf)
    if len(ag) < len(pf):
        agi = np.arange(len(ag))
        np.random.shuffle(agi)
        pfi = np.random.choice(np.arange(len(pf)), len(ag), replace=False)
    else:
        pfi = np.arange(len(pf))
        np.random.shuffle(pfi)
        agi = np.random.choice(np.arange(len(ag)), len(pf), replace=False)

    # data for pairs
    u1, v1, u2, v2 = get_vectors(ag, pf, agi, pfi)
    a1, a2 = get_angles(ag, pf, agi, pfi)
    w1, w2 = get_mags(ag, pf, agi, pfi)

    # test correlation
    vc = veccor(u1, v1, u2, v2)
    ac = stats.spearmanr(a1,a2)[0]
    wc = stats.spearmanr(w1,w2)[0]

    vec_correlation[P] = vc
    ang_correlation[P] = ac
    mag_correlation[P] = wc

# shove into dataframe and save
permutation = pd.DataFrame({'vector':vec_correlation, 'angle':ang_correlation, 'mag':mag_correlation})
permutation.to_csv(args.out+'_permutation.txt', sep='\t')
