import pandas as pd, numpy as np, os, glob, re

def mean_err(plocs, metadata):
    plocs = pd.read_csv(plocs)
    plocs = plocs.rename(columns={'x':'xpred','y':'ypred'})
    metadata = pd.read_csv(metadata, sep='\t')
    plocs = pd.merge(plocs, metadata, on='sampleID')

    xerr = np.mean(plocs.xpred - plocs.x)
    yerr = np.mean(plocs.ypred - plocs.y)
    merr = np.mean(np.hypot((plocs.xpred - plocs.x),
                            (plocs.ypred - plocs.y)))

    return xerr, yerr, merr

predlocs = glob.glob('out/predlocs/*predlocs.txt')
x_error = np.empty(len(predlocs))
y_error = np.empty(len(predlocs))
m_error = np.empty(len(predlocs))
sigmas = np.empty(len(predlocs))
skews = np.empty(len(predlocs))

for i, pred in enumerate(predlocs):
    if 'uniform' in pred:
        sigma, mu, run = re.findall(r'\d+(?:\.\d+)?', pred)
        skew = 0
    else:
        sigma, mu, run, skew = re.findall(r'\d+(?:\.\d+)?', pred)
    meta = f'out/metadata/sigma_{sigma}_bias_{mu}_run_{run}_metadata.txt'
    x, y, m = mean_err(pred, meta)
    
    x_error[i] = x
    y_error[i] = y
    m_error[i] = m
    sigmas[i] = sigma
    skews[i] = skew

df = pd.DataFrame({'x_err':x_error,
                   'y_err':y_error,
                   'm_err':m_error,
                   'sigma':sigmas,
                   'skew':skews})
df.to_csv('simulation_error_450train_50test_unweighted.txt', sep='\t')
