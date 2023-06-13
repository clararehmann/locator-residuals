# spatially downsample training set using GRTS geosampling
# https://github.com/jgrss/geosample

import argparse
import pandas as pd, numpy as np, geopandas as gpd
from geosample import QuadTree

parser=argparse.ArgumentParser()
parser.add_argument('--data', help='path to metadata TSV with "sampleID", "x", and "y" columns')
parser.add_argument('--outpath', help='path to training set output TSV (appended with "_GRTSsample.txt")')
parser.add_argument('--train_size', default=10000, help='maximum number of individuals to keep (will be much lower, default 10000)')
args=parser.parse_args()

def GRTS_sample(train_size, test_ind, metadata):
    # make GeoDataFrame
    geometa = gpd.GeoDataFrame(metadata, geometry=gpd.points_from_xy(metadata.x, metadata.y))
    # GRTS sample
    qt = QuadTree(geometa)
    qt.split_recursive(first_null=True)
    d = qt.sample(train_size)
    loc_df = d[['sampleID','x','y']].copy()
    # test set
    test = metadata.loc[test_ind].copy()
    # dataframe of train and test individuals
    loc_df = pd.concat([loc_df, test]).sample(frac=1).reset_index(drop=True)
    return loc_df

# read input metadata, identify individuals to use in test set
metadata = pd.read_csv(args.data, sep='\t', index_col=0)
test_ind = np.isnan(metadata.x)
loc_df = GRTS_sample(args.train_size, test_ind, metadata)
loc_df.to_csv(args.outpath+'_GRTSsample.txt', sep='\t')
