import pandas as pd

sra_df = pd.read_csv(snakemake.input[1])
sra_df = sra_df[["Run", "GEO_Accession"]]
# I needed to manually remove "!" from the begining of several lines
df = pd.read_csv(snakemake.input[0], sep="\t", comment="!",
                 header=None, index_col=0)
df = df.transpose()
# filter for only RNA-seq samples
df = df.loc[df["Sample_description"] == 'single cell RNA seq']
# remove HeLa contaminants - maybe use as a positive control?
df = df.loc[~(df['Sample_characteristics_ch2'] == "tissue: HeLa contaminant")]
cols_to_keep = ["Sample_title", "Sample_geo_accession"]
df = df[cols_to_keep]
df = df.merge(sra_df, left_on="Sample_geo_accession", right_on="GEO_Accession")
df = df.rename(columns={"Run": "sample"})
# next step - get patient information, site, etc. from Sample_title
df["patient_id"] = df["Sample_title"].apply(lambda x: x.split("_")[1])
df["location"] = df["Sample_title"].apply(lambda x: x.split("_")[2])
df["patient"] = df["Sample_title"].apply(lambda x: "{}_{}".format(x.split("_")[1], x.split("_")[2]))
patients = df[["patient", "patient_id", "location"]]
patients.drop_duplicates().to_csv(snakemake.output[0], index=False, sep="\t")
samples = df.drop(columns=["Sample_geo_accession", "patient_id", "location"])
samples.drop_duplicates().to_csv(snakemake.output[1], index=False, sep="\t")
# df.to_csv("data/sample-metadata.tsv", sep="\t", index=False)
