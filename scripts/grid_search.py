import pandas as pd, numpy as np
import sys, os, glob, argparse, itertools

parser=argparse.ArgumentParser()
parser.add_argument('--lossdir', help='path to directory of training histories formatted TS_lambda_LM_bandwidth_BD_history.txt')
parser.add_argument('--out', help='path to save tsv of final training losses (appended with _grid_search.txt')
args = parser.parse_args()

bnds = [10, 100, 1000, 10000]
lamb = ['0.1', '0.5', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

def get_loss(b, l):
    loss = []
    for i in range(10):
        if os.path.exists(f'{args.lossdir}/{i}_lambda_{l}_bandwidth_{b}_history.txt'):
            df = pd.read_csv(f'{args.lossdir}/{i}_lambda_{l}_bandwidth_{b}_history.txt', sep='\t')
            ls = np.array(df['val_loss'])[-1]
            loss.append(ls)
    if len(loss) > 0:
        return np.mean(loss)
    else:
        return None

grid = itertools.product(*[bnds, lamb])
grid = np.array(list(grid))

validation_loss = []

for b, l in grid:
    validation_loss.append(get_loss(b, l))

df = pd.DataFrame({'bandwidth': grid[:,0],
                    'lambda': grid[:,1],
                    'val_loss': validation_loss})
df.to_csv(args.out + '_grid_search.txt', sep='\t')
