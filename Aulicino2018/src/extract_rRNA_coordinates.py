import pandas as pd



col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

df = pd.read_csv(snakemake.input[0], sep="\t", names=col_names, comment="#")

if snakemake.wildcards["filter"] == "rRNA":
# let's extract all rRNA reads
    rRNA_df = df.loc[df.attribute.str.contains("gene_biotype \"rRNA\"")]
    rRNA_df[["seqname", "start", "end"]].to_csv(snakemake.output[0], sep="\t", index=False, header=False)

if snakemake.wildcards["filter"] == "16S":
    rRNA_16S_df = df.loc[df.attribute.str.contains("product \"16S ribosomal RNA\"")]
    rRNA_16S_df[["seqname", "start", "end"]].to_csv(snakemake.output[0], sep="\t", index=False, header=False)

if snakemake.wildcards["filter"] == "protein_coding":
# let's extract all rRNA reads
    protein_coding_df = df.loc[df.attribute.str.contains("gene_biotype \"protein_coding\"")]
    protein_coding_df[["seqname", "start", "end"]].to_csv(snakemake.output[0], sep="\t", index=False, header=False)
