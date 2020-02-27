import pandas as pd

# we care about time, status, infection
df = pd.read_csv(snakemake.input[0])

patients = df[['time', 'status', 'infection']].drop_duplicates()
patients = patients.reset_index(drop=True)
patients["patient"] = ["Pt{}".format(i) for i in patients.index]

samples = df[['time', 'status', 'infection', 'Run', 'plate', 'well']]
samples = samples.rename(columns={'Run': 'sample'})
samples = samples.merge(patients, on=['time', 'status', 'infection'])

patients.to_csv(snakemake.output[0], sep="\t", index=False)
samples.drop(columns=['time', 'status', 'infection']).to_csv(snakemake.output[1], sep="\t", index=False)
