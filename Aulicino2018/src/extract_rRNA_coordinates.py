import pandas as pd



col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

df = pd.read_csv("data/GCF_000006945.2_ASM694v2_genomic.gtf", sep="\t", names=col_names, comment="#")

if snakemake.wildcards["filter"] == "rRNA":
# let's extract all rRNA reads
    rRNA_df = df.loc[df.attribute.str.contains("gene_biotype \"rRNA\"")]
    rRNA_df.to_csv(snakemake.output[0], columns=["seqname", "start", "stop"], sep="\t", index=False, header=False)

if snakemake.wildcards["filter"] == "16S":
    rRNA_16S_df = df.loc[df.attribute.str.contains("product \"16S ribosomal RNA\"")]
    rRNA_16S_df.to_csv(snakemake.output[0], columns=["seqname", "start", "stop"], sep="\t", index=False, header=False)
