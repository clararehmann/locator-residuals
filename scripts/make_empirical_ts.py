# make training sets from an empirical dataset,
# use every individual for prediction

import argparse
import numpy as np, pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--metadata', help='path to sample data tsv with "sampleID", "x", and "y" columns')
parser.add_argument('--pred_prop', type=float, help='proportion of individuals to hold out as a prediction set')
parser.add_argument('--out', help='outpath for training datasets')
args=parser.parse_args()

meta=pd.read_csv(args.metadata,sep='\t')
ind=np.arange(len(meta))
np.random.shuffle(ind)
na_count=int(len(meta)/((1-args.pred_prop)*len(meta)))
groups=np.array_split(ind,na_count)
for g in range(len(groups)):
    meta=pd.read_csv(args.metadata,sep='\t')
    for i in groups[g]:
        meta.loc[i,'x']=np.nan
        meta.loc[i,'y']=np.nan
    meta.to_csv(args.out+'_'+str(g)+'.txt',sep='\t')
