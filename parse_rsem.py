import os
import pandas as pd

sample_df = pd.read_csv('metadata_all_studies.csv')
samples = sample_df['sample_id']
files = [os.path.join('rsem', f"{sample}-{condition}.isoforms.results") for sample, condition in zip(sample_df.sample_id, sample_df.condition)]
# files = [os.path.join('rsem', sample, '.isoforms.results') for sample in samples]
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

