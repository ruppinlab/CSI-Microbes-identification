import pandas as pd
import pysam
from collections import defaultdict

# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile(snakemake.input[0], mode="rb")

df = pd.read_csv(snakemake.input[2], sep="\t", header=None)
barcodes = df[0].values

d = defaultdict(list)
# seg is an AlignedSegment object
for seg in pathseq_bam.fetch(until_eof=True):
    cb = seg.query_name.split("_")[1]
    d[cb].append(seg)

# we want to write an output bam for every key in the dictionary
# cellbarcodes = d.keys()
output_bam_file = "output/PathSeq/{patient}-{sample}-{plate}-{cell}/pathseq.bam"
for cellbarcode in barcodes:
    f = output_bam_file.format(patient=snakemake.wildcards["patient"], sample=snakemake.wildcards["sample"], plate=snakemake.wildcards["plate"], cell=cellbarcode)
    output_bam = pysam.AlignmentFile(f, mode="wb", template=pathseq_bam)
    for seg in d[cellbarcode]:
        output_bam.write(seg)
    output_bam.close()

pathseq_bam.close()
