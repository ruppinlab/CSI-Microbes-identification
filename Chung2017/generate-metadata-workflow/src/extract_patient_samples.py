import pandas as pd

df = pd.read_csv(snakemake.input[0])
df = df.rename(columns={"Run": "sample", "patient_id": "patient"})
patient_cols = ['Cell_type', 'patient', 'molecular_subtype']
patients = df[patient_cols]
sample_cols = ["sample", "patient"]
samples = df[sample_cols]
patients.to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
