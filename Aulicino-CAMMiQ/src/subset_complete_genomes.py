import pandas as pd

column_names = [
    "assembly_accession", "bioproject", "biosample", "wgs_master",
    "refseq_category", "taxid", "species_taxid", "organism_name",
    "infraspecific_name", "isolate", "version_status", "assembly_level",
    "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter",
    "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq",
    "relation_to_type_material"
    ]

df = pd.read_csv(snakemake.input[0], sep="\t", skiprows=2, names=column_names)

# only use complete genomes
df = df.loc[df["assembly_level"] == "Complete Genome"]

# in the future, exclude the genomes that are excluded from RefSeq (only 4/13737)
# df = df.loc[df["excluded_from_refseq"].isna()]

# only keep the first genome for each species
df = df.drop_duplicates(subset="species_taxid")

# write the entire file
df.to_csv(snakemake.output[0], sep="\t", index=False)
# write the ftp_path's to file so we can use for rsync command
