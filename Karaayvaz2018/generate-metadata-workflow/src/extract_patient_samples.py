import pandas as pd

df = pd.read_csv(snakemake.input[0])
patient_cols = ["patient", "breast_cancer_subtype", "tissue"]
patients = df[patient_cols]
sample_cols = ["Run", "patient", "Bases"]
samples = df[sample_cols]
samples = samples.rename(columns={"Run": "sample"})
patients.to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
