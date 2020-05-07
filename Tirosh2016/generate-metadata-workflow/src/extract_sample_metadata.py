import pandas as pd
import re
import numpy as np

# cell_df = pd.read_csv("data/celltype_annotation.tsv", sep="\t")
# cell_df = cell_df.iloc[1:]
# cell_df = cell_df.rename(columns={"NAME": "sample", "CLUSTER": "celltype", "SUB-CLUSTER": "celltype1"})
# cell_df["Tumor"] = cell_df["celltype"].apply(lambda x: "Tumor" if x == "malignant" else "nonTumor")

cell_df = pd.read_csv("data/GSE72056_melanoma_single_cell_revised_v2.txt", sep="\t", nrows=3, index_col=0).transpose()
cell_df = cell_df.rename(columns={
    "tumor": "patient",
    "malignant(1=no,2=yes,0=unresolved)": "malignant",
    "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)": "celltype"
    })
cell_df["malignant"] = cell_df["malignant"].replace({1: "nonmalignant", 2: "malignant", 0: "unknown"})
cell_df["celltype"] = cell_df.apply(lambda x: "malignant" if x.malignant == "malignant" else x.celltype, axis=1)
cell_df["celltype"] = cell_df["celltype"].replace({
    1: "Tcell",
    2: "Bcell",
    3: "Macrophage",
    4: "Endothelial",
    5: "CAF",
    6: "NK",
    0: "unknown"
    })
cell_df["patient"] = cell_df["patient"].apply(lambda x: "Melanoma_{}".format(x))
cell_df = cell_df.reset_index()
cell_df = cell_df.rename(columns={"index": "sample"})
sample_df = pd.read_csv("data/sample.tsv", sep="\t")
# drop participant because it is frequently wrong
# ex. cy72-CD45-pos-H12-S960 belongs to Melanoma 74 - one sample
# ex. all 224 cy94-* assigned to Melanoma 75 but all CY94 samples are assigned to Melanoma 94
sample_df = sample_df.drop(columns="participant")
sample_df = sample_df.rename(columns={'entity:sample_id': 'sample'})

patient_df = pd.read_csv("data/participant.tsv", sep="\t")
patient_df = patient_df.rename(columns={"entity:participant_id": "patient"})
# of the 4645 cells from 19 patient, only 1396 from 5 patients have overlapping sample names
# some (Melanoma_65, Melanoma_67) - has the weird sample names like 2_HA5P2ADXX_2_AACATAAT_ACACGATC
# some (Melanoma_67) start with Cy67 and end with _L001_R1_001_fastq_gz
# many end with -comb or _comb, which differs across the cell and sample dataframes
# this increases to 1484 from 6 patients
sample_df["sample"] = sample_df["sample"].str.replace("_L001_R1_001_fastq_gz", "")
# ~2000 cell_df["sample"] end with -comb
# ~2400 cell_df["sample"] end with _comb
# 0 sample_df["sample"] end with "-comb"
# 1396 sample_df["sample"] end with "_comb"
# if we replace "-comb" and "_comb" with "", we still get 4645 unique sample names from cell_df and sample_df
# so we aren't losing the uniqueness of any samples
cell_df["sample"] = cell_df["sample"].str.replace("_comb", "").str.replace("-comb", "")
sample_df["sample"] = sample_df["sample"].str.replace("_comb", "").str.replace("-comb", "")
# now when we merge the samples, we get 4505 cells from 17 patients
# Melanoma_59 - cell_df["sample"] starts with tumor1-cell but sample_df["sample"] starts with Cy59
# sample_df["sample"] ends with _45_S45 but cell_df["sample"] ends with _45
pat = r"_(?P<sample>\d{1,2})"
repl = lambda m: "_S{}".format(m.group('sample'))
sample_df.loc[sample_df["sample"].str.startswith("tumor1"), "sample"] = sample_df.loc[sample_df["sample"].str.startswith("tumor1-cell-"), "sample"].str.replace("tumor1-cell-", "Cy59_").str.replace("_S(\d)+", "").str.replace(pat, repl)
cell_df.loc[cell_df["sample"].str.startswith("Cy59"), "sample"] = cell_df.loc[cell_df["sample"].str.startswith("Cy59"), "sample"].str.replace(pat, repl)
# now we get 4575 cells from 18 patients - we are content with this for now
sample_df = sample_df.merge(cell_df, on="sample")
# how many sample names end with something like S291?
# all sample names do now
sample_df["sample_num"] = sample_df["sample"].str.extract("S(\d{1,4})")
sample_df["sample_batch"] = np.floor(sample_df.sample_num.astype("int").subtract(1).divide(96)).astype("int")
# how many sample names end with something like A03-S291?
# sometimes, the sample-derived patient and the listed patient don't agree

sample_df["sample"] = sample_df["sample"].str.replace("-", "_")
# if we try to split to extract the "patient", it doesn't work too well
# df["sample"].apply(lambda x: x.split("_")[0]).unique()
# returns values like 'cy88cd45', 'cy88', 'cy89a', 'cy89core11', 'cy89coreq1', etc.
sample_df.to_csv(snakemake.output[1], sep="\t", index=False)
patient_df.loc[patient_df.patient.isin(sample_df.patient)].to_csv(snakemake.output[0], sep="\t", index=False)
