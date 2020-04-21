import pandas as pd
import re
# there are 10 patients
# 9/10 patients have a mixture of tumor and non-tumor cells
# patient MGH103 contains only tumor cells
# for each patient with both types of cells,
# we have one file containing all the tumor samples
# and one file containing all the non-tumor samples
# samples are of the form {patient}_{plate}_{well}
# let's ignore MGH103 for now
patients = [
    "MGH42", "MGH43", "MGH44", "MGH45", "MGH56", "MGH57", "MGH61", "MGH64",
    "MGH103", "MGH107"
    ]

output_df = []

for patient in patients:
    tumor_samples = pd.read_csv("data/{}_tumor_columns_selected.txt".format(patient), sep="\t", nrows=0, index_col=0).columns
    tumor_df = pd.DataFrame(data={"patient": patient, "sample": tumor_samples, "celltype": "Tumor"})
    output_df.append(tumor_df)
    if patient != "MGH103":
        other_samples = pd.read_csv("data/{}_other_columns_selected.txt".format(patient), sep="\t", nrows=0, index_col=0).columns
        other_df = pd.DataFrame(data={"patient": patient, "sample": other_samples, "celltype": "nonTumor"})
        output_df.append(other_df)

df = pd.concat(output_df)
df["sample"] = df["sample"].apply(lambda x: x if x.startswith("MGH") else "MGH{}".format(x))
df["sample"] = df["sample"].str.replace("_", "-")
df["plate"] = df["sample"].apply(lambda x: x.split("-")[1]) #x.split("_")[1])
df["well"] = df["sample"].apply(lambda x: x.split("-")[2])
df.to_csv(snakemake.output[0], sep="\t", index=False)
df["cancertype"] = "Astrocytoma"
df[["patient", "cancertype"]].drop_duplicates().to_csv(snakemake.output[1], sep="\t", index=False)
