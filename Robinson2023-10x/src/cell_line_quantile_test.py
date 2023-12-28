import pandas as pd

# kmer_df = pd.read_csv("output/SAHMI/P1/SCAF2963_3_Live.barcode.kmer.hits.tsv", sep="\t")
kmer_df = pd.read_csv(snakemake.input["kmer"], sep="\t")
# report_df = pd.read_csv("output/SAHMI/P1/SCAF2963_3_Live.kraken.report.rpmm.tsv", sep="\t")
report_df = pd.read_csv(snakemake.input["rpmm"], sep="\t")
# cell_lines_df = pd.read_excel("../../SAHMI/Table S4.xlsx")
cell_lines_df = pd.read_excel(snakemake.input["cell_line"])
merged_df = kmer_df.merge(report_df.drop(columns=["study", "sample"]), on=["taxid", "name"])
merged_df = merged_df.loc[merged_df.p.notna()]
cell_lines_df = cell_lines_df[["name", "rank", "taxid", 0.99]]
merged_df = merged_df.merge(cell_lines_df, on=["name", "rank", "taxid"])
hits_df = merged_df.loc[merged_df.rpmm > merged_df[0.99]]
hits_df.to_csv(snakemake.output["hits"], sep="\t")
hits_df["taxid"].to_csv(snakemake.output["taxa"], sep="\t", header=False, index=False)
