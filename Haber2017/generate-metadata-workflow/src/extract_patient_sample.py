import pandas as pd

def get_infection(x):
    if isinstance(x.treatment, str):
        if "Salmonella" in x.treatment:
            return "Salmonella"
        if "polygyrus" in x.treatment:
            return "H. polygyrus"
    if isinstance(x.agent, str):
        if "Salmonella" in x.agent:
            return "Salmonella"
    return "None"

df = pd.read_csv("data/SraRunTable.csv")
df["seq_type"] = ["10x" if "TenX" in x else "Smart-seq2" for x in df["DATASTORE filetype"]]

# treatment and agent should be one column
df["infection"] = df.apply(get_infection, axis=1)

patients = df[["Organism", "source_name", "infection"]].drop_duplicates()
patients = patients.reset_index(drop=True)
patients["patient"] = ["Pt{}".format(i) for i in patients.index]
patients.to_csv(snakemake.output[0], sep="\t", index=False)

samples = df[["Run", "Instrument", "seq_type", "Bases", "Instrument", "LibraryLayout", "Organism", "source_name", "infection"]]
samples = samples.rename(columns={'Run': 'sample'})
samples = samples.merge(patients, on=["Organism", "source_name", "infection"])
samples = samples.drop(columns=["Organism", "source_name", "infection"])
samples.to_csv(snakemake.output[1], sep="\t", index=False)
