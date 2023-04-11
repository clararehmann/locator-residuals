import pandas as pd, numpy as np
import os, subprocess

"""
download AG1000G v3 data
"""

# download metadata
#os.makedirs('v3/metadata')
#subprocess.run('gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ v3/metadata', shell=True)

# download zarrs
os.makedirs('v3/snp_genotypes/all')
subprocess.run('gsutil -m rsync -r gs://vo_agam_release/v3/snp_genotypes/all v3/snp_genotypes/all', shell=True)

# sort metadata - only gambiae
general_path = 'v3/metadata/general/'
species_path = 'v3/metadata/species_calls_aim_20220528/'

cohorts = os.listdir(species_path)

ag1000g_metadata = pd.DataFrame()

for cohort in cohorts:
    data = os.path.join(general_path, cohort)
    metadata = pd.read_csv(os.path.join(data, 'samples.meta.csv'))
    snpdata = pd.read_csv(os.path.join(data, 'wgs_snp_data.csv'))
    data = pd.merge(metadata, snpdata, on=list(set(metadata.columns) & set(snpdata.columns)))
    species = pd.read_csv(os.path.join(species_path, cohort, 'samples.species_aim.csv'))
    data = pd.merge(data, species, on=list(set(data.columns) & set(species.columns)))
    data = data[data['aim_species_gambiae_coluzzii'] == 'gambiae']
    ag1000g_metadata = pd.concat([ag1000g_metadata, data])

ag1000g_metadata.to_csv('v3/ag1000g_all_samples_metadata.txt', sep='\t')
