import pandas as pd
import numpy as np


def return_celltype(x):
    if x == "unknown":
        return "unknown"
    if x == "epithelial":
        return "Tumor"
    return "nonTumor"


sra_df = pd.read_csv(snakemake.input[0])
sra_df = sra_df.rename(columns={"Run": "sample", "GEO_Accession (exp)": "GSM"})
celltype_df = pd.read_csv(snakemake.input[1], sep="\t")
celltype_df.columns = ["sample_name", "celltype"]
# use gse_df to convert between SRA sample and celltype sample
gse_df = pd.read_csv(snakemake.input[2], sep="\t", header=None).transpose()
gse_df.columns = ["sample_name", "GSM"]
# use this file to extract "batch_number" and "batch_depleted"
qc_df = pd.read_csv(snakemake.input[3], sep="\t")
qc_df = qc_df[["number_batch", "depletion_batch"]]
post_QC_df = pd.read_csv(snakemake.input[4], sep="\t")
# merge all the files together
df = celltype_df.merge(gse_df, on="sample_name", how="right")
df = df.merge(sra_df, on="GSM", how="right")
df = df.merge(qc_df, left_on="sample_name", right_index=True)
patient_cols = ["patient", "breast_cancer_subtype", "tissue"]
patients = df[patient_cols]
sample_cols = ["sample", "patient", "sample_name", "celltype", "number_batch", "depletion_batch"]
samples = df[sample_cols]
samples = samples.rename(columns={"Run": "sample", "number_batch": "batch"})
samples["celltype"] = samples["celltype"].fillna("unknown")
samples["cancer"] = samples["celltype"].apply(return_celltype)
samples["QC_status"] = samples["sample_name"].apply(lambda x: "Passed" if x in post_QC_df.Sample else "Failed")
# samples["plate"] = samples["sample_name"].apply(lambda x: x.split("_")[1])
samples["plate"] = samples.apply(lambda x: "{}-{}".format(x["batch"], x["sample_name"].split("_")[1]), axis=1)
#print(samples.plate)
samples["well"] = samples["sample_name"].apply(lambda x: x.split("_")[2])
patients.drop_duplicates().to_csv(snakemake.output[0], index=False, sep="\t")
samples.to_csv(snakemake.output[1], index=False, sep="\t")
