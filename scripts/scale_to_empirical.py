import pandas as pd, numpy as np
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--centroids', help='path to SLiM centroid data to rescale')
parser.add_argument('--scale', help='path to x, y scale factor data (--scout from sample_like_empirical.py)')
parser.add_argument('--slimlandscape', default=50, help='maximum coordinates of SLiM landscape')
parser.add_argument('--out', help='path to write rescaled centroid data to')
args=parser.parse_args()

def rescale(old_value, old_min, old_max, new_min, new_max):
    new_value = ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    return new_value

# read in dataframes
slimdf = pd.read_csv(args.centroids, sep='\t')
scaledf = pd.read_csv(args.scale, sep='\t', index_col=0)

# get scale factors
xmin, xmax = scaledf.loc['min','x'], scaledf.loc['max','x']
ymin, ymax = scaledf.loc['min','y'], scaledf.loc['max','y']

# rescale SLiM data to empirical coordinates
reslim = slimdf.copy()

for i in range(len(reslim)):
    x=reslim.loc[i, 'x']
    y=reslim.loc[i, 'y']
    gc_x=reslim.loc[i, 'gc_x']
    gc_y=reslim.loc[i, 'gc_y']
    kd_x=reslim.loc[i, 'kd_x']
    kd_y=reslim.loc[i, 'kd_y']

    reslim.loc[i, 'x']  = rescale(x, 0, args.slimlandscape, xmin, xmax)
    reslim.loc[i, 'y'] = rescale(y, 0, args.slimlandscape, ymin, ymax)
    reslim.loc[i, 'gc_x'] = rescale(gc_x, 0, args.slimlandscape, xmin, xmax)
    reslim.loc[i, 'gc_y'] = rescale(gc_y, 0, args.slimlandscape, ymin, ymax)
    reslim.loc[i, 'kd_x'] = rescale(kd_x, 0, args.slimlandscape, xmin, xmax)
    reslim.loc[i, 'kd_y'] = rescale(kd_y, 0, args.slimlandscape, ymin, ymax)

reslim.to_csv(args.out+'.txt', sep='\t')
