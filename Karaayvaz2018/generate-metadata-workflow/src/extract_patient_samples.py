import pandas as pd
import numpy as np


def return_celltype(x):
    if x == "QC-filtered":
        return "QC-filtered"
    if x == "epithelial":
        return "Tumor"
    return "nonTumor"


sra_df = pd.read_csv(snakemake.input[0])
sra_df = sra_df.rename(columns={"Run": "sample", "GEO_Accession (exp)": "GSM"})
celltype_df = pd.read_csv(snakemake.input[1], sep="\t")
celltype_df.columns = ["sample_name", "celltype"]
gse_df = pd.read_csv(snakemake.input[2], sep="\t", header=None).transpose()
gse_df.columns = ["sample_name", "GSM"]
df = celltype_df.merge(gse_df, on="sample_name", how="right")
df = df.merge(sra_df, on="GSM", how="right")
patient_cols = ["patient", "breast_cancer_subtype", "tissue"]
patients = df[patient_cols]
sample_cols = ["sample", "patient", "sample_name", "celltype"]
samples = df[sample_cols]
samples = samples.rename(columns={"Run": "sample"})
samples["celltype"] = samples["celltype"].fillna("QC-filtered")
samples["cancer"] = samples["celltype"].apply(return_celltype)
samples["plate"] = samples["sample_name"].apply(lambda x: x.split("_")[1])
patients.drop_duplicates().to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
