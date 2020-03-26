import pandas as pd

sra_df = pd.read_csv(snakemake.input[0])
sra_df = sra_df.rename(columns={"Run": "sample", "GEO_Accession (exp)": "GSM"})
celltype_df = pd.read_csv(snakemake.input[1], sep="\t")
celltype_df.columns = ["sample_name", "celltype"]
gse_df = pd.read_csv(snakemake.input[2], sep="\t", header=None).transpose()
gse_df.columns = ["sample_name", "GSM"]
df = celltype_df.merge(gse_df, on="sample_name")
df = df.merge(sra_df, on="GSM")
patient_cols = ["patient", "breast_cancer_subtype", "tissue"]
patients = df[patient_cols]
sample_cols = ["sample", "patient", "sample_name", "celltype"]
samples = df[sample_cols]
samples = samples.rename(columns={"Run": "sample"})
patients.to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
