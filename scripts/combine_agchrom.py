import allel, re, os, matplotlib, sys, zarr, time, subprocess, copy
import numpy as np, pandas as pd
import dask

nsites = 1000000

# script to put all Anopheles genome data in one Zarr file
# for grid search

def process_chromosome(gt, samples, nsites):
    gindex = np.argsort(samples)

    sites = np.sort(np.random.choice(range(len(gt)), nsites, replace=False))
    genotypes = gt[sites,:,:]
    genotypes = genotypes[:,gindex,:]

    genotypes = allel.GenotypeArray(genotypes)

    tmp = genotypes.count_alleles()
    biallel=tmp.is_biallelic()

    genotypes=genotypes[biallel,:,:]

    derived_counts=genotypes.count_alleles()[:,1]
    ac_filter=[x >= 2 for x in derived_counts]
    genotypes=allel.GenotypeDaskArray(genotypes[ac_filter,:,:])
    
    return genotypes

chromosome = '2L.zarr'
gt = zarr.open('snp_genotypes_all/' + chromosome + '/calldata/GT')
samples= zarr.open('snp_genotypes_all/' + chromosome + '/samples')
#positions = zarr.open('snp_genotypes_all/' + chromosome + '/variants/POS' )
gt = allel.GenotypeDaskArray(gt)

genotypes = process_chromosome(gt, samples, nsites)
print(chromosome + ' complete')

for chromosome in ['2R.zarr','3L.zarr','3R.zarr','X.zarr']:
    gt_ = zarr.open('snp_genotypes_all/' + chromosome + '/calldata/GT')
    gt_ = allel.GenotypeDaskArray(gt_)
    samples_ = zarr.open('snp_genotypes_all/'+chromosome+'/samples')
    genotypes_ = process_chromosome(gt_, samples_, nsites)
    genotypes = dask.array.concatenate([genotypes, genotypes_], axis=0)
    print(chromosome + ' complete')

genotypes = genotypes[np.sort(np.random.choice(range(len(genotypes)), 100000, replace=False)), :, :]    
genotypes = genotypes.rechunk()
samples = np.sort(samples)

store = zarr.DirectoryStore('snp_genotypes_all/ALL_CHROM')
root = zarr.group(store=store)
calldata = root.create_group('calldata')
#samples = root.create_group('samples')
dask.array.to_zarr(genotypes,'snp_genotypes_all/ALL_CHROM/calldata/GT', overwrite=True)
root.create_dataset(name='samples', data=samples, shape=len(samples))
#dask.array.to_zarr(samples, 'snp_genotypes_all/ALL_CHROM/samples', overwrite=True)
#dask.array.to_zarr(positions, 'snp_genotypes_all/ALL_CHROM/variants/POS', overwrite=True)    
