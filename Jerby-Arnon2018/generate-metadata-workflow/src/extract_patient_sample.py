import pandas as pd

patient_df = pd.read_csv(snakemake.input[0], sep="\t", skiprows=1)
patient_df = patient_df.rename(columns={"default_participant": "patient"})
sample_df = pd.read_csv(snakemake.input[1], sep="\t")
sample_df = sample_df.rename(columns={"participant": "patient", "entity:sample_id": "sample"})

sample_df["plate"] = sample_df["sample"].apply(lambda x: x.split("_")[0])

patient_df.to_csv(snakemake.output[0], sep="\t", index=False)
sample_df.to_csv(snakemake.output[1], sep="\t", index=False)
