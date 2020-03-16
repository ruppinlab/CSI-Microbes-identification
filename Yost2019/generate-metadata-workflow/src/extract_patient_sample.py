import pandas as pd

df = pd.read_csv(snakemake.input[0])
patients = df[["Donor", "treatment", "Cell_type"]].drop_duplicates()
patients = patients.reset_index(drop=True)
patients["patient"] = ["Pt{}".format(i) for i in patients.index]

samples = df[["Run", "Assay Type", "AvgSpotLen", "Donor", "treatment", "Cell_type"]]
samples = samples.rename(columns={'Run': 'sample'})
samples = samples.merge(patients, on=["Donor", "treatment", "Cell_type"])

patients.to_csv(snakemake.output[0], sep="\t", index=False)
samples.drop(columns=["Donor", "treatment", "Cell_type"]).to_csv(snakemake.output[1], sep="\t", index=False)
