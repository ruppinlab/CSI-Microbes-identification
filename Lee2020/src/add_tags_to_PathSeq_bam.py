import pysam
import pandas as pd
from collections import defaultdict

# load and index the CellRanger BAM file
cr_bam = pysam.AlignmentFile(pathseq.input[0], mode="rb")
cr_idx = pysam.IndexedReads(cr_bam)
cr_idx.build()

# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile(pathseq.input[1], mode="rb")

iter = pathseq_bam.fetch()
output = []
d = defaultdict(lambda: defaultdict(list))
# seg is an AlignedSegment object
for seg in iter:
    # returns an IteratorRowSelection object, which contains one or more AlignedSegment object
    cr_list = list(cr_idx.find(seg.query_name))
    # we assume that all records belonging to the same query name will have the same CB/UB tag
    # not all records will have the CB tag and the UB tag
    if cr_list[0].has_tag("CB") and cr_list[0].has_tag("UB"):
        CB = cr_list[0].get_tag(tag="CB")
        UB = cr_list[0].get_tag(tag="UB")
        # using set_tags removes all other tags - use set_tag instead
        seg.set_tag("CB", CB, "Z")
        seg.set_tag("UB", UB, "Z")
        d[CB][UB].append(seg)
    # keep all PathSeq alignments
    output.append(seg)

# write all PathSeq alignments with or without tags
all_pathseq_bam = pysam.AlignmentFile(snakemake.output[0], mode="wb", template=pathseq_bam)
for seg in output:
    all_pathseq_bam.write(seg)

all_pathseq_bam.close()

# desired output is one BAM file per valid cell barcode with only one read per UMI
# read in the units dataframe
# subset by the sample wildcards
df = pd.read_csv(snakemake.input[2], sep="\t")
df = df.loc[(df["patient"] == snakemake.wildcards["patient"]) and (df["sample"] == snakemake.wildcards["sample"])]
#barcodes = ["AAGGAGCGTGAGTATA-1", "ACATGGTAGCAGATCG-1", "ACGATACGTCTCTCGT-1", "ACGTCAAGTGTTGAGG-1"]
barcodes = df["barcode"]
# cell_iter = d["AAGGAGCGTGAGTATA-1"]
out_file = "output/PathSeq/{}-{}-{}/pathseq_with_tags.bam".format(snakemake.wildcards["patient"], snakemake.wildcards["sample"], "{}")
for barcode in barcodes:
    UMI_dict = d[barcode]
    barcode_bam = pysam.AlignmentFile(out_file.format(barcode), mode="wb", template=pathseq_bam)
    for UMI in UMI_dict:
        # keep one read per UMI - the read with the highest mapping quality
        UMI_reads = UMI_dict[UMI]
        UMI_read = UMI_reads[0]
        for read in UMI_reads:
            if read.mapping_quality > UMI_read.mapping_quality:
                UMI_read = read
        barcode_bam.write(UMI_read)
    barcode_bam.close()

cr_bam.close()
pathseq_bam.close()
