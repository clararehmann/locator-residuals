import pandas as pd, numpy as np
import zarr, os, dask, allel, subprocess

agpath = '/sietch_colab/crehmann/vo_agam_release/v3'

# sort metadata - only gambiae mosquitoes
cohorts = os.listdir(os.path.join(agpath, 'metadata/species_calls_aim_20220528/'))
ag1000g_metadata = pd.DataFrame()
for cohort in cohorts:
    metadata = pd.read_csv(os.path.join(agpath, 'metadata/general', cohort, 'samples.meta.csv'))
    snpdata = pd.read_csv(os.path.join(agpath, 'metadata/general', cohort, 'wgs_snp_data.csv'))
    spcdata = pd.read_csv(os.path.join(agpath, 'metadata/species_calls_aim_20220528', cohort, 'samples.species_aim.csv'))
    data = pd.merge(metadata, snpdata, on=['sample_id'])
    data = pd.merge(data, spcdata, on=['sample_id'])
    data = data[data['aim_species_gambiae_coluzzii'] == 'gambiae']

    ag1000g_metadata = pd.concat([ag1000g_metadata, data])

ag1000g_metadata = pd.DataFrame({'sampleID':ag1000g_metadata.sample_id,
                                 'x':ag1000g_metadata.longitude,
                                 'y':ag1000g_metadata.latitude})

ag1000g_metadata.to_csv('data/ag1000g/ag1000g_v3_gambiae.txt', sep='\t')


# combine genotype data across cohorts
if os.path.exists(os.path.join(agpath, 'snp_genotypes_all')):
    pass
else:
    os.makedirs(os.path.join(agpath, 'snp_genotypes_all'))

def concat_zarr(gt, samples, cohort, chromosome):
    gt_ = zarr.open(os.path.join(agpath, 'snp_genotypes/all', cohort, chromosome, 'calldata/GT'))
    samples_ = zarr.open(os.path.join(agpath, 'snp_genotypes/all', cohort, 'samples'))
    gt_ = allel.GenotypeDaskArray(gt_)
    gt = dask.array.concatenate([gt, gt_], axis=1)
    samples = dask.array.concatenate([samples, samples_])
    return gt, samples

for chromosome in ['2L','2R','3L','3R','X']:
    # copy positions over
    os.makedirs(os.path.join(agpath, 'snp_genotypes_all', chromosome, 'variants'), exist_ok=True)
    ps = f'cp -r {os.path.join(agpath, "snp_genotypes/all/sites", chromosome, "variants/POS")} {os.path.join(agpath, "snp_genotypes_all", chromosome, "variants")}'
    subprocess.run(ps, shell=True)
 
    print(f'concatenating chromosome {chromosome}')
    cohorts = os.listdir(os.path.join(agpath, 'metadata/species_calls_aim_20220528'))
    gt = zarr.open(os.path.join(agpath, 'snp_genotypes/all', cohorts[0], chromosome, 'calldata/GT'))
    gt = allel.GenotypeDaskArray(gt)
    samples = zarr.open(os.path.join(agpath, 'snp_genotypes/all', cohorts[0], 'samples'))
    print(f'{cohorts[0]} complete')
    cohorts.pop(0)
    while len(cohorts) > 0:
        gt, samples = concat_zarr(gt, samples, cohorts[0], chromosome)
        print(f'{cohorts[0]} complete')
        cohorts.pop(0)
    gt = gt.rechunk()
    samples = samples.rechunk()    
    store = zarr.DirectoryStore(os.path.join(agpath, 'snp_genotypes_all', chromosome))
    root = zarr.group(store=store)
    calldata = root.create_group('calldata')
    dask.array.to_zarr(gt, os.path.join(agpath, 'snp_genotypes_all', chromosome, 'calldata/GT'))
    dask.array.to_zarr(samples, os.path.join(agpath, 'snp_genotypes_all', chromosome, 'samples'))

