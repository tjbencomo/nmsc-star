import os
import pandas as pd


samples = os.listdir('rsem')
files = [os.path.join('rsem', sample, '.isoforms.results') for sample in samples]
# files = [f for f in files if os.path.isfile(f)]

dfList = []
for sample, f in zip(samples, files):
    print(f"Loading {sample}")
    df = pd.read_csv(f, sep = '\t')
    df = df[['transcript_id', 'TPM']]
    df['SampleID'] = sample
    dfList.append(df)
df = pd.concat(dfList, ignore_index = True)
print(df.head())
df.to_csv('rsem_tpm.txt.gz', index = False)

