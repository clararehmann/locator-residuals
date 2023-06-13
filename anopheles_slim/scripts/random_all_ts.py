import argparse
import numpy as np, pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--metadata')
parser.add_argument('--ts',type=float,default=0.9)
parser.add_argument('--use_all',default=False,action="store_true",help="cycle thru and make multiple ts files")
parser.add_argument('--out')
args=parser.parse_args()

meta=pd.read_csv(args.metadata,sep='\t')

if args.use_all:
    ind=np.arange(len(meta))
    np.random.shuffle(ind)
    na_count=int(len(meta)/((1-args.ts)*len(meta)))
    groups=np.array_split(ind,na_count)
    for g in range(len(groups)):
        meta=pd.read_csv(args.metadata,sep='\t')
        for i in groups[g]:
            meta.loc[i,'x']=np.nan
            meta.loc[i,'y']=np.nan
        meta.to_csv(args.out+'_'+str(g)+'.txt',sep='\t')
else:
    drop=np.random.choice(range(len(meta)),int((1-args.ts)*len(meta)),replace=False) # randomly choose samples to NaN
    for i in drop:
        meta.loc[i,'x']=np.nan
        meta.loc[i,'y']=np.nan
    meta.to_csv(args.out,sep='\t')
