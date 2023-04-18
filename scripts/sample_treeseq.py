# process tree sequence,
# save full metadata for all individuals,
# create training sets (uniform, biased, point-sampled)

from slimtools import *
import subprocess

training_sets = ['uniform', 'point', 'half']
Nreps = 1
skews = [0.5, 0.6, 0.7, 0.8, 0.9]
tssize = 450

outpath = 'out/'
metadatapath = 'out/metadata/'
vcfpath = 'out/genotypes/'
tspath = 'out/training_sets/'
for f in [outpath, metadatapath, vcfpath, tspath]:
    subprocess.run(f'mkdir -p {f}', shell=True)

tree = sys.argv[1]
sim, sigma, mu = np.array(tree.split('/')[-1].replace('.trees','').split('_'))[[0,2,4]].astype(float)
sim = int(sim)

# process each tree sequence 10x
for n in range(Nreps):
    ts = tskit.load(tree)
    ts = process_treeseq(ts)
    
    all_samples = []
    
    # sample uniformly
    uniform_df = sample_ts(ts, 'uniform', tssize, max_distance=1.5)
    all_samples.extend(uniform_df.sampleID)
    uniform_df.to_csv(os.path.join(tspath, f'{sim}_sigma_{sigma}_mu_{mu}_uniform_{n}.txt'), sep='\t')
    
    for skew in skews:
        # point sampling
        point_df = sample_ts(ts, 'point', tssize, point=(40,40), decay=1/skew)
        all_samples.extend(point_df.sampleID)
        point_df.to_csv(os.path.join(tspath, f'{sim}_sigma_{sigma}_mu_{mu}_point_{skew}_{n}.txt'), sep='\t')
        
        # biased sampling
        half_df = sample_ts(ts, 'half', tssize, skew=skew)
        all_samples.extend(half_df.sampleID)
        half_df.to_csv(os.path.join(tspath, f'{sim}_sigma_{sigma}_mu_{mu}_half_{skew}_{n}.txt'), sep='\t')

    # simplify tree sequence to only relevant individuals
    all_samples = np.unique(all_samples)
    individual_ids = [int(i.replace('tsk_','')) for i in all_samples]
    nodes = []
    for i in individual_ids:
        nodes.extend(ts.individual(i).nodes)
    ts = ts.simplify(nodes)
    
    # save individual locations
    metadata = sample_locations(ts)
    metadata.to_csv(os.path.join(metadatapath, f'{sim}_sigma_{sigma}_mu_{mu}_{n}.txt'), sep='\t')
    
    # save genotypes
    treepath = os.path.join(vcfpath, f'{sim}_sigma_{sigma}_mu_{mu}_{n}.vcf')
    
    with open(treepath, 'w') as f:
        ts.write_vcf(f, individual_names=all_samples)
