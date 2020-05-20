import pandas as pd


df = pd.read_csv("data/SraRunTable.txt")
df = df.rename(columns={"Run": "sample"})
df["patient"] = "Pt0"
df = df[["sample", "patient", "infection"]]
df[["patient"]].drop_duplicates().to_csv("output/patients.tsv", sep="\t", index=False)
df.to_csv("output/samples.tsv", sep="\t", index=False)
