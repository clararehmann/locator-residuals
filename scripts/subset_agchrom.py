import allel, re, os, matplotlib, sys, zarr, time, subprocess, copy
import numpy as np, pandas as pd
import dask

nsites = 200000
chromosome = '2R'

# subset N SNPs from a chromosome,
# save as Zarr

def process_chromosome(gt, nsites):
    sites = np.sort(np.random.choice(range(len(gt)), nsites, replace=False))
    genotypes = gt[sites,:,:]

    genotypes = allel.GenotypeArray(genotypes)

    tmp = genotypes.count_alleles()
    biallel=tmp.is_biallelic()

    genotypes=genotypes[biallel,:,:]

    derived_counts=genotypes.count_alleles()[:,1]
    ac_filter=[x >= 2 for x in derived_counts]
    genotypes=allel.GenotypeDaskArray(genotypes[ac_filter,:,:])
    
    return genotypes

gt = zarr.open('snp_genotypes_all/' + chromosome + '.zarr/calldata/GT')
samples= zarr.open('snp_genotypes_all/' + chromosome + '.zarr/samples')
positions = zarr.open('snp_genotypes_all/' + chromosome + '.zarr/variants/POS' )
gt = allel.GenotypeDaskArray(gt)

genotypes = process_chromosome(gt, nsites*20)
while len(genotypes) < nsites:
    gt2 = process_chromosome(gt, nsites*20)
    genotypes = dask.array.concatenate([genotypes, gt2], axis=0)
sites = np.sort(np.random.choice(range(len(genotypes)), nsites, replace=False))
genotypes = genotypes[sites,:,:]
genotypes = genotypes.rechunk()

store = zarr.DirectoryStore('snp_genotypes_all/'+chromosome+'_subset_'+str(nsites))
root = zarr.group(store=store)
calldata = root.create_group('calldata')
variants = root.create_group('variants')
#samples = root.create_group('samples')
dask.array.to_zarr(genotypes,'snp_genotypes_all/'+chromosome+'_subset_'+str(nsites)+'/calldata/GT', overwrite=True)
root.create_dataset(name='samples', data=samples, shape=len(samples))
variants.create_dataset(name='POS', data=positions, shape=positions.shape)
#dask.array.to_zarr(samples, 'snp_genotypes_all/ALL_CHROM/samples', overwrite=True)
#dask.array.to_zarr(positions, 'snp_genotypes_all/ALL_CHROM/variants/POS', overwrite=True)    
