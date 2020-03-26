import pandas as pd

sra_df = pd.read_csv(snakemake.input[0])
sra_df = sra_df.rename(columns={"Run": "sample", "patient_id": "patient", "GEO_Accession (exp)": "GSM"})
gse_df = pd.read_csv(snakemake.input[1], sep="\t", header=None).transpose()
gse_df.columns = ["GSM", "sample_name"]
celltype_df = pd.read_csv(snakemake.input[2], sep="\t")
celltype_df = celltype_df.loc[celltype_df["type"] == "SC"]
celltype_df = celltype_df.rename(columns={"sample": "sample_name", "index": "celltype1", "index2": "celltype2", "index3": "celltype3"})
df = gse_df.merge(celltype_df, on="sample_name")
df = df.merge(sra_df, on="GSM")
patient_cols = ['Cell_type', 'patient', 'molecular_subtype']
patients = df[patient_cols].drop_duplicates()
sample_cols = ["sample", "patient", "sample_name", "celltype1", "celltype2", "celltype3"]
samples = df[sample_cols]

patients.to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
