import pandas as pd

sra_df = pd.read_csv(snakemake.input[0])
sra_df = sra_df.rename(columns={"Run": "sample", "patient_id": "patient", "GEO_Accession (exp)": "GSM", 'Cell_type': 'site'})
gse_df = pd.read_csv(snakemake.input[1], sep="\t", header=None).transpose()
gse_df.columns = ["GSM", "sample_name"]
celltype_df = pd.read_csv(snakemake.input[2], sep="\t")
celltype_df = celltype_df.loc[celltype_df["type"] == "SC"]
celltype_df = celltype_df.rename(columns={"sample": "sample_name", "index": "Tumor", "index2": "celltype1", "index3": "celltype2"})
df = gse_df.merge(celltype_df, on="sample_name")
df = df.merge(sra_df, on="GSM")
df["plate"] = df["patient"]
df["patient"] = df.patient.str.replace("LN", "")
df["patient"] = df.patient.str.replace("_Re", "")
patient_cols = ['patient', 'molecular_subtype']
patients = df[patient_cols].drop_duplicates()
sample_cols = ["sample", "patient", "sample_name", "site",  "Tumor", "celltype1", "celltype2", "plate"]
samples = df[sample_cols]

patients.to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
