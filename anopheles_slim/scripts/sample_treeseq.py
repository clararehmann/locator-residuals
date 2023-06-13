## load in mutated (or raw) tree sequence, sample Anopheles analogs,
## write metadata and zarr

import msprime, pyslim, argparse, allel, os, sys
from scipy.spatial import KDTree
from itertools import chain
import numpy as np, pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--tree', help='path to .trees file from SLiM output')
parser.add_argument('--process', default=False, action='store_true', help='recapitate and add neutral mutations?')
parser.add_argument('--out', help='output stem (will be appended with .zarr and _metadata.txt')
parser.add_argument('--scalefactor', type=float, help='factor by which SLiM africa landscape is scaled')
parser.add_argument('--anopheles', default='/home/crehmann/kernlab/Locator/anopheles/ag1000g_metadata.txt', help='path to Anopheles metadata')
parser.add_argument('--rrate', type=float, default=1e-8, help='recombination rate')
parser.add_argument('--mrate', type=float, default=1e-8, help='mutation rate')
args=parser.parse_args()

def rescale(old_value, old_min, old_max, new_min, new_max):
    new_value = ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    return new_value

def make_KDtree(x, y):
    xyval = np.stack((x, y), axis=1)
    tree=KDTree(xyval)
    return tree

# load and process data
ts=pyslim.load(args.tree)

if args.process:
    rts=ts.recapitate(recombination_rate=args.rrate) # recapitate at same recombination rate
    last_gen=rts.individuals_alive_at(0)
    last_nodes=[]
    for i in last_gen:
        last_nodes.extend(rts.individual(i).nodes)
    sts=rts.simplify(last_nodes)
    ts=pyslim.SlimTreeSequence(msprime.mutate(sts,rate=args.mrate,keep=True)) # add neutral mutations

# get SLiMulation metadata (unique SLiM ID, x and y coordinates)    
locs=ts.individual_locations[:, :2]
xval = locs[:, 0]
yval = locs[:, 1]
slimid = [i.metadata['pedigree_id'] for i in ts.individuals()]

# set up scaling
### Anopheles
ag = pd.read_csv(args.anopheles, sep='\t')
axmin, axmax = -16.5, 44.7
aymin, aymax = -27, 15.4
ag = ag[[axmin <= x <= axmax for x in ag.x]]
### SLiM
sxmin, sxmax = 0, 6900*args.scalefactor
symin, symax = 0, 4748*args.scalefactor

# rescale SLiMulation data for KDTree sampling
re_x = [rescale(x, sxmin, sxmax, axmin, axmax) for x in xval]
re_y = [rescale(y, symin, symax, aymin, aymax) for y in yval]

# set up SLiMulation metadata
agslim = pd.DataFrame({'x':re_x, 'y':re_y, 'slim':slimid, 'pair_x':np.repeat(np.nan, len(re_x)), 'pair_y':np.repeat(np.nan, len(re_y))})
agslim['tree0_index'] = agslim.index  # save original tskit index
workingdf = agslim.copy() # remove sampled SLiM individuals as we go (no repeats when querying the tree again)
agslim.index = agslim.slim # set SLiM ID as index in original dataframe

# find SLiMulation nearest neighbors of each Anopheles datapoint
emprxy = np.stack((ag.x, ag.y), axis=1)
tree = make_KDtree(agslim.x, agslim.y)

for coord in np.random.permutation(np.unique(emprxy, axis=0)): # for unique anopheles sampling location
    
    # get same number of SLiM samples near that location
    s = len(ag[(ag.x == coord[0]) & (ag.y == coord[1])])
    q = tree.query(coord, k=s)[1]
    
    # add paired Anopheles point to original SLiM dataframe
    if type(q) == np.int64:
        slimID = workingdf.loc[q, 'slim']
        agslim.loc[slimID, 'pair_x'] = coord[0]
        agslim.loc[slimID, 'pair_y'] = coord[1]
    else:
        for v in q:
            slimID = workingdf.loc[v, 'slim']
            agslim.loc[slimID, 'pair_x'] = coord[0]
            agslim.loc[slimID, 'pair_y'] = coord[1]
    
    # remove samples from working dataframe so they're not used again
    workingdf.drop(q, inplace=True)
    workingdf.reset_index(inplace=True, drop=True)
    tree = make_KDtree(workingdf.x, workingdf.y)

# individuals to keep in the tree
keep_ids = agslim.slim[~np.isnan(agslim.pair_x)]
keep_inds = list(agslim.tree0_index[~np.isnan(agslim.pair_x)])

# simplify tree
keep_nodes = []
for i in keep_inds:
    keep_nodes.extend(ts.individual(i).nodes)
ts = ts.simplify(keep_nodes)

# deal with 'missing' data (https://tskit.dev/pyslim/docs/latest/tutorial.html#extracting-individuals-after-simplification)
indivlist = []
for i in ts.individuals_alive_at(0):
    ind = ts.individual(i)
    if ts.node(ind.nodes[0]).is_sample():
        indivlist.append(i)
        # if one node is a sample, the other should be also:
        assert ts.node(ind.nodes[1]).is_sample()

# metadata for new tree
sampled_agslim = agslim.loc[keep_ids].reset_index(drop=True) # cut down dataframe to saved individuals
sampled_agslim.index = sampled_agslim.slim # make SLiM ID the index
sampled_agslim['tree1_index'] = np.repeat(np.nan, len(keep_ids))
# save simplified tree sequence indices
for individual in indivlist:
    sampled_agslim.loc[ts.individual(individual).metadata['pedigree_id'], 'tree1_index'] = int(individual)
# make sure everything is in order!
assert list(sampled_agslim.tree1_index) == sorted(sampled_agslim.tree1_index)

# final formatted metadata for Locator
metadata = pd.DataFrame({'sampleID':[str(s) for s in sampled_agslim.slim], 'x':sampled_agslim.pair_x, 'y':sampled_agslim.pair_y}).reset_index(drop=True)

# write vcf
with open(args.out+'.vcf','w') as vcf:
     ts.write_vcf(vcf, individuals=indivlist, individual_names=list(metadata.sampleID))

# write metadata
metadata.to_csv(args.out+'_metadata.txt', sep='\t')

# convert to zarr
allel.vcf_to_zarr(args.out+'.vcf',args.out+'.zarr')
