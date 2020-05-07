import pandas as pd

df = pd.read_csv("data/E-MTAB-5061.sdrf.txt", sep="\t")

df = df.rename(columns={
    'Source Name' : "sample",
    "Characteristics[single cell well quality]": "well_quality",
    "Characteristics[cell type]": "cell_type",
    "Comment[FASTQ_URI]": "fastq_url",
    "Characteristics[individual]": "patient",
    "Characteristics[sex]": "sex",
    "Characteristics[disease]": "disease"
})

sample_df = df[["sample", "well_quality", "cell_type", "fastq_url", "patient"]]
patient_df = df[["patient", "sex", "disease"]].drop_duplicates()

sample_df.to_csv("output/samples.tsv", sep="\t", index=False)
patient_df.to_csv("output/patients.tsv", sep="\t", index=False)

# we want
# for sample_df
# sample - 'Source Name'
# well quality - 'Characteristics[single cell well quality]'
# cell_type - 'Characteristics[cell type]'
# fastq_url - "Comment[FASTQ_URI]"
# for patient_df
# patient - 'Characteristics[individual]'
# sex - 'Characteristics[sex]'
# disease - 'Characteristics[disease]'
