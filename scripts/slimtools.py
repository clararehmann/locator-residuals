import sys, msprime, numpy as np, pandas as pd
import tskit, pyslim, os, io, shutil, itertools
import scipy, allel
from scipy import spatial

def process_treeseq(infile,
                   recombination_rate=1e-8,
                   mutation_rate=1e-8,
                   seed=0):
    """
    recapitates, mutates, and simplifies tree sequence to the last generation of individuals
    returns tree sequence and table of sample locations
    """
    if type(infile) == str:
        ts = pyslim.load(infile)
    else:
        ts = infile
    if(not seed==0):
        np.random.seed(seed)
    # recapitate
    ts = pyslim.recapitate(ts, 
                           recombination_rate=recombination_rate,
                           ancestral_Ne=len(pyslim.individuals_alive_at(ts, 0)))
    # simplify to last generation
    last_gen = pyslim.individuals_alive_at(ts, 0)
    last_nodes = []
    for i in last_gen:
        last_nodes.extend(ts.individual(i).nodes)
    ts = ts.simplify(last_nodes)
    # add neutral variation
    ts = msprime.sim_mutations(ts,
                               rate=mutation_rate,
                               model=msprime.SLiMMutationModel(type=0),
                               keep=True)
    return ts

def sample_locations(ts):
    """
    get locations of all samples in the tree sequence
    """
    locs = ts.tables.individuals.location
    locs = np.reshape(locs, (len(ts.tables.individuals), 3))
    locs = locs[:,0:2]
    samples = [f'tsk_{i}' for i in np.arange(len(ts.tables.individuals))]
    df = pd.DataFrame({'sampleID':samples,
                       'x':locs[:, 0],
                       'y':locs[:,1]})
    return df

def uniform_sample(ts, N, max_distance, exclude_ind=None):
    # get samples and locations
    locs = ts.tables.individuals.location
    locs = np.reshape(locs, (len(ts.tables.individuals), 3))
    locs = locs[:,0:2]
    samples = np.array([f'tsk_{i}' for i in np.arange(len(ts.tables.individuals))])
    if exclude_ind is not None:
        locs = np.delete(locs, exclude_ind, axis=0)
        samples = np.delete(samples, exclude_ind)
    
    # grid to uniformly sample over
    x = np.linspace(0, 50, int(np.sqrt(N+1)))
    y = np.linspace(0, 50, int(np.sqrt(N+1)))
    xy = np.array(list(itertools.product(x, y)))
    
    # randomly sample a proximate individual to each node
    dists = spatial.distance.cdist(xy, locs)
    sample_inds = np.argsort(dists)
    dist_mask = np.sort(dists) < max_distance
    sample_inds = [sample_inds[i][dist_mask[i]] for i in range(len(sample_inds))]
    sample_inds = np.array([np.random.choice(i) if len(i) > 0 else np.nan for i in sample_inds])
    sample_inds = sample_inds[~np.isnan(sample_inds)]
    sample_inds = sample_inds.astype(int)
    sample_inds = np.unique(sample_inds.astype(int))
    
    samples = samples[sample_inds]
    locs = locs[sample_inds]
    
    return sample_inds, samples, locs

def point_sample(ts, N, point, decay, exclude_ind=None):
    locs = ts.tables.individuals.location
    locs = np.reshape(locs, (len(ts.tables.individuals), 3))
    locs = locs[:,0:2]
    samples = np.array([f'tsk_{i}' for i in np.arange(len(ts.tables.individuals))])
    if exclude_ind is not None:
        locs = np.delete(locs, exclude_ind, axis=0)
        samples = np.delete(samples, exclude_ind)
    
    dists_to_point = [scipy.spatial.distance.euclidean(locs[x],point) for x in range(len(locs))]
    weights = [1/x**decay for x in dists_to_point]
    weights = weights/np.sum(weights)
    
    sample_inds = np.random.choice(range(len(weights)), N, p=weights, replace=False)
    samples = samples[sample_inds]
    locs = locs[sample_inds]
    
    return samples, locs

def half_sample(ts, N, skew, midpoint=25, exclude_ind=None):
    locs = ts.tables.individuals.location
    locs = np.reshape(locs, (len(ts.tables.individuals), 3))
    locs = locs[:,0:2]
    samples = np.array([f'tsk_{i}' for i in np.arange(len(ts.tables.individuals))])
    if exclude_ind is not None:
        locs = np.delete(locs, exclude_ind, axis=0)
        samples = np.delete(samples, exclude_ind)
    
    weights = np.empty(len(samples))
    weights[np.where(locs[:,0] > midpoint)] = skew
    weights[np.where(locs[:,0] <= midpoint)] = 1-skew
    weights = weights/np.sum(weights)
    
    sample_inds = np.random.choice(range(len(weights)), N, p=weights, replace=False)
    samples = samples[sample_inds]
    locs = locs[sample_inds]
    
    return samples, locs

def sample_ts(ts,
              sampling,
              training_set,
              prediction_set=50,
              point=(45,45),
              skew=0.5,
              decay=0.5,
              max_distance=1
              ):
    """
    sample a uniformly-distributed prediction set and 
    [uniformly, point, or half-skewed] training set from a tree sequence
    ts: tree sequence object
    sampling: sampling method for training set
    training_set: number of individuals to include in 
    sampling methods: 'uniform': ,'point','half'
    point: location of target point for point sampling
    bias: degree of bias for sampling
    """
    if sampling not in ['uniform','point','half']:
        sys.exit()
    else:
        print(sampling, training_set)
        # uniformly sample prediction set
        pred_ind, pred_samp, pred_locs = uniform_sample(ts, prediction_set, max_distance)
        print(len(pred_ind))
        # sample training set
        if sampling=='uniform':
            train_inds, train_samp, train_locs = uniform_sample(ts, training_set, max_distance, exclude_ind=pred_ind)
        elif sampling=='point':
            train_samp, train_locs = point_sample(ts, training_set, point, 1/skew, exclude_ind=pred_ind)
        elif sampling=='half':
            train_samp, train_locs = half_sample(ts, training_set, skew, midpoint=25, exclude_ind=pred_ind)
    
    # combine two, make dataframe
    samples = np.concatenate([pred_samp, train_samp])
    locs = np.concatenate([np.full(pred_locs.shape, np.nan), train_locs])
    df = pd.DataFrame({'sampleID':samples, 
                       'x':locs[:,0],
                       'y':locs[:,1]})
    
    # simplify tree sequence
#    individual_ids = [int(i.replace('tsk_','')) for i in samples]
#    nodes = []
#    for i in individual_ids:
#        nodes.extend(ts.individual(i).nodes)
#    ts = ts.simplify(nodes)
    
    return df#, ts
