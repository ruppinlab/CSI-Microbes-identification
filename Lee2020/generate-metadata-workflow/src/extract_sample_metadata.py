import pandas as pd

df = pd.read_csv("data/E-MTAB-8410.sdrf.txt", sep="\t")

column_name_table = {
    "Comment[FASTQ_URI]": "cell_barcode_fastq_url",
    "Comment[FASTQ_URI].1": "cDNA_fastq_url",
    "Comment[FASTQ_URI].2": "IDX_fastq_url",
    "Characteristics[sampling site]": "sample_site",
    'Source Name': "sample",
    "Characteristics[sex]": "sex",
    "Characteristics[individual]": "patient",
    "Characteristics[age]": "age"
}

df = df.rename(columns=column_name_table)
df["lane"] = df["Assay Name"].apply(lambda x: x.split("_")[3])
patients = df[["patient", "age", "sex"]].drop_duplicates()
samples = df[["sample", "sample_site", "lane", "patient",
              "cell_barcode_fastq_url", "cDNA_fastq_url", "IDX_fastq_url"]]

patients.to_csv("output/patients.tsv", sep="\t", index=False)
samples.to_csv("output/samples.tsv", sep="\t", index=False)
