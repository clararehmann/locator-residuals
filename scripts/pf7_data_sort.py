import pandas as pd, numpy as np, os, subprocess

# only individuals that pass QC and are in africa
xlim = (-31.3, 58)
ylim = (-37.5, 38.8)

pfpath = '/sietch_colab/crehmann/pf7/Pf7_samples.txt'
outpath = '.'

df = pd.read_csv(pfpath, sep='\t')
df = df[df['QC pass'] == True].reset_index(drop=True)
df['x'] = df['Admin level 1 longitude']
df['y'] = df['Admin level 1 latitude']
df = df[(xlim[0] < df.x) & (xlim[1] > df.x) & (ylim[0] < df.y) & (ylim[1] > df.y)]
df = pd.DataFrame({'sampleID':df.Sample,
                    'x':df.x,
                    'y':df.y})
df.to_csv(os.path.join(outpath, 'pf7_africa_QC.txt'), sep='\t')

