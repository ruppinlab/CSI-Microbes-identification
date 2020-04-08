import pandas as pd

# we care about time, status, infection
df = pd.read_csv(snakemake.input[0])

patients = df[['time', 'status', 'infection']].drop_duplicates()
patients = patients.reset_index(drop=True)
# patients["old_patient"] = ["Pt{}".format(i) for i in patients.index]

samples = df[['Run', 'time', 'status', 'infection', 'plate', 'well']]
samples = samples.rename(columns={'Run': 'sample'})
samples = samples.loc[samples["status"].isin(["exposed", "infected", "uninfected"])]
samples = samples.merge(patients, on=['time', 'status', 'infection'])
samples["patient"] = "Pt0"
patients["patient"] = "Pt0"
samples["infected"] = samples["status"].apply(lambda x: "uninfected" if x == "uninfected" else "infected")
patients.drop(columns=['time', 'status', 'infection']).drop_duplicates().to_csv(snakemake.output[0], sep="\t", index=False)
samples.to_csv(snakemake.output[1], sep="\t", index=False)
