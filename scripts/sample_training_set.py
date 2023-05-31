import numpy as np, pandas as pd
import argparse, os, io, scipy, itertools, zarr
from scipy import spatial

parser=argparse.ArgumentParser()
parser.add_argument('--metadata')
parser.add_argument('--zarr')
parser.add_argument('--ts_size', default=450)
parser.add_argument('--pred_size', default=50)
parser.add_argument('--max_distance', default=2)
parser.add_argument('--out')
args = parser.parse_args()

df = pd.read_csv(args.metadata, sep='\t')
zarr = zarr.open_group(args.zarr)
samples = zarr['samples'][:]
df = df[[i in samples for i in df['sampleID']]]

def uniform_sample(df, N, max_distance, exclude_ind=None):
    locs = df[['x','y']].to_numpy()
    samples = df.sampleID.to_numpy()
    if exclude_ind is not None:
        locs = np.delete(locs, exclude_ind, axis=0)
        samples = np.delete(samples, exclude_ind)
    
    # grid to uniformly sample over
    x = np.linspace(0, 50, int(np.sqrt(N+1)))
    y = np.linspace(0, 50, int(np.sqrt(N+1)))
    xy = np.array(list(itertools.product(x,y)))
    
    # randomly sample a proximate individual to each node
    dists = spatial.distance.cdist(xy, locs)
    sample_inds = np.argsort(dists)
    dist_mask = np.sort(dists) < max_distance
    sample_inds = [sample_inds[i][dist_mask[i]] for i in range(len(sample_inds))]
    sample_inds = np.array([np.random.choice(i) if len(i) > 0 else np.nan for i in sample_inds])
    sample_inds = sample_inds[~np.isnan(sample_inds)]
    sample_inds = sample_inds.astype(int)
    sample_inds = np.unique(sample_inds.astype(int))
    
    return locs[sample_inds], samples[sample_inds], sample_inds
    
def half_sample(df, N, skew, midpoint=25, exclude_ind=None):
    locs = df[['x','y']].to_numpy()
    samples = df.sampleID.to_numpy()
    if exclude_ind is not None:
        locs = np.delete(locs, exclude_ind, axis=0)
        samples = np.delete(samples, exclude_ind)
    
    weights = np.empty(len(samples))
    weights[np.where(locs[:,0] > midpoint)] = skew
    weights[np.where(locs[:,0] <= midpoint)] = 1-skew
    weights = weights/np.sum(weights)

    sample_inds = np.random.choice(range(len(weights)), N, p=weights, replace=False)
    return locs[sample_inds], samples[sample_inds], sample_inds

# get prediction samples
plocs, psamp, pinds = uniform_sample(df, args.pred_size, args.max_distance)
pred = pd.DataFrame({'sampleID': psamp,
                     'x': np.full(len(psamp), np.nan),
                     'y': np.full(len(psamp), np.nan)})

# sample uniform training set
locs, samps, inds = uniform_sample(df, args.ts_size, args.max_distance, pinds)
ts = pd.DataFrame({'sampleID': samps,
                   'x': locs[:,0],
                   'y': locs[:,1]})
ts = pd.concat((pred, ts))
ts = ts.sample(frac = 1)
ts.to_csv(args.out+'_uniform_ts.txt', sep='\t')

for skew in [0.6, 0.7, 0.8, 0.9]:
    locs, samps, inds = half_sample(df, args.ts_size, skew, midpoint=25, exclude_ind=pinds)
    ts = pd.DataFrame({'sampleID': samps,
                       'x': locs[:,0],
                       'y': locs[:,1]})
    ts = pd.concat((pred, ts))
    ts = ts.sample(frac = 1)
    ts.to_csv(args.out +'_skew_' + str(skew) + '_ts.txt', sep='\t')
