import pandas as pd

df = pd.read_csv("data/SraRunTable.txt")

df = df.rename(columns={"Run": "sample", "sequencing_batch": "plate"})
df["patient"] = "Pt0"
sample_df = df[["sample", "Condition", "Hours_post_infection", "plate", "patient"]]
patient_df = df[["patient"]]
sample_df.to_csv("output/samples.tsv", sep="\t", index=False)
patient_df.drop_duplicates().to_csv("output/patients.tsv", sep="\t", index=False)
