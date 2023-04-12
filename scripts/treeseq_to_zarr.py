import msprime, pyslim, argparse, allel, os
import numpy as np, pandas as pd
from matplotlib import pyplot as plt

parser=argparse.ArgumentParser(description='write last generation of SLiM output to zarr, save location metadata for each individual')
parser.add_argument('--tree', help='path to .trees file from SLiM output')
parser.add_argument('--out', help='output stem (will be appended with .zarr and _metadata.txt')
parser.add_argument('--rrate', type=float, default=1e-8, help='recombination rate')
parser.add_argument('--mrate', type=float, default=1e-8, help='mutation rate')
parser.add_argument('--plot', default=False, action='store_true', help='plot locations of all individuals on landscape')
parser.add_argument('--output', type=str, default='zarr', choices=['zarr', 'vcf'], help='type of output genotype file (default zarr, options vcf or zarr)')
args=parser.parse_args()

#def lastgen(treeseq): # simplify tree sequence to last generation
#    last_gen=treeseq.individuals_alive_at(0)
#    keep_inds=np.random.choice(last_gen,args.n_size,replace=False)
#    last_nodes=[]
#    for i in keep_inds:
#        last_nodes.extend(treeseq.individual(i).nodes)
#    treeseq=treeseq.simplify(last_nodes)
#    return treeseq

#def write_file(treeseq): # write vcf and metadata
#    with open(args.out+'.vcf','w') as vcf:
#        treeseq.write_vcf(vcf)
#    sampleID=[]
#    x=[]
#    y=[]
#    for ind in treeseq.individuals():
#        sampleID.append('tsk_'+str(ind.id))
#        x.append(ind.location[0])
#        y.append(ind.location[1])
#    meta=pd.DataFrame(list(zip(sampleID,x,y)),columns=['sampleID','x','y'])
#    meta.to_csv(args.out+'_metadata.txt',sep='\t')
#    return meta

def plot(df):
    fig,ax=plt.subplots()
    ax.scatter(df.x,df.y,alpha=0.1,s=25)
    ax.set_aspect('equal')
    plt.savefig(args.out+'.png')
    return None

# load and process data
ts=pyslim.load(args.tree)
rts=ts.recapitate(recombination_rate=args.rrate) # recapitate at same recombination rate
last_gen=rts.individuals_alive_at(0)
last_nodes=[]
for i in last_gen:
    last_nodes.extend(rts.individual(i).nodes)
sts=rts.simplify(last_nodes)
fts=pyslim.SlimTreeSequence(msprime.mutate(sts,rate=args.mrate,keep=True)) # add neutral mutations

# write vcf
with open(args.out+'.vcf','w') as vcf:
    fts.write_vcf(vcf, individuals=fts.individuals_alive_at(0))
sampleID=[]
x=[]
y=[]
for ind in fts.individuals():
    sampleID.append('tsk_'+str(ind.id))
    x.append(ind.location[0])
    y.append(ind.location[1])
meta=pd.DataFrame(list(zip(sampleID,x,y)),columns=['sampleID','x','y'])
meta.to_csv(args.out+'_metadata.txt',sep='\t')



# write data

#meta=write_file(fts) # vcf and metadata
#if args.plot: # plot if you want
#    plot(meta)
if args.output == 'zarr': # write to zarr if you want
    allel.vcf_to_zarr(args.out+'.vcf',args.out+'.zarr') 
