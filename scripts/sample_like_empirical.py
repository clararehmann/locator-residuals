# sample NxN spatial SLiMulation with same spatial distribution as empirical dataset

import pandas as pd, numpy as np
from scipy.spatial import KDTree
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--slim', help='path to SLiM (x,y) metadata for selective sampling')
parser.add_argument('--slimlandscape', default=50, help='maximum coordinates of SLiM landscape')
parser.add_argument('--mimic', help='path to empirical (x,y) metadata to mimic during sampling')
parser.add_argument('--scout', help='path to write x, y scaling factors to (appended with _scalefactors.txt). (necessary to scale back to mimic data for comparison). columns = x, y      rows = min, max')
parser.add_argument('--lssize', default=None, type=float, help='NxN square landscape size of empirical data for SLiM scaling. default is built on scale of data')
parser.add_argument('--out', help='path to write selectively-sampled SLiM data to')
args=parser.parse_args()

def rescale(old_value, old_min, old_max, new_min, new_max):
    new_value = ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    return new_value

# read in dataframes
slimdf = pd.read_csv(args.slim, sep='\t')
emprdf = pd.read_csv(args.mimic, sep='\t')

# find center of empirical data
xsize = max(emprdf.x) - min(emprdf.x) 
xcenter = max(emprdf.x) - (xsize/2)
ysize = max(emprdf.y) - min(emprdf.y)
ycenter =  max(emprdf.y) - (ysize/2)

# rescale empirical data to SLiM coordinates
empslim = emprdf.copy()

if args.lssize:
    lssize=args.lssize
else:
    lssize=max(xsize,ysize) + (args.slimlandscape/5)
xmin = xcenter - (lssize/2)
xmax = xcenter + (lssize/2)
ymin = ycenter - (lssize/2)
ymax = ycenter + (lssize/2)

scaledf = pd.DataFrame([[xmin,ymin],[xmax,ymax]], columns=['x','y'], index=['min','max'])
scaledf.to_csv(args.scout+'_scalefactors.txt',sep='\t')

for i in range(len(empslim)):
    x=empslim.loc[i, 'x']
    y=empslim.loc[i, 'y']

    empslim.loc[i, 'x']  = rescale(x, xmin, xmax, 0, args.slimlandscape)
    empslim.loc[i, 'y'] = rescale(y, ymin, ymax, 0, args.slimlandscape)

# sample SLiM data by nearest neighbors

slimxy = np.stack((slimdf.x, slimdf.y), axis=1)
emprxy = np.stack((empslim.x, empslim.y), axis=1)
tree = KDTree(slimxy)

inds = []
xs = []
ys = []

for coord in np.unique(emprxy, axis=0):
    s = len(empslim[(empslim.x == coord[0]) & (empslim.y == coord[1])])
    q = tree.query(coord, s)[1]
    if type(q) == np.int64: 
        inds.append(q)
        xs.append(coord[0])
        ys.append(coord[1])
    else:
        inds.extend(q)
        for i in range(len(q)):
            xs.append(coord[0])
            ys.append(coord[1])

# save individuals to dataframe

mimic = pd.DataFrame(columns = ['sampleID', 'x','y'])
for i in range(len(inds)):
    x = xs[i]
    y = ys[i]
    ID = slimdf.iloc[inds[i]].sampleID
    mimic = pd.concat([mimic, pd.DataFrame({'sampleID':[ID], 'x':[x], 'y':[y]})], ignore_index=True)

mimic.to_csv(args.out+'.txt', sep='\t')
