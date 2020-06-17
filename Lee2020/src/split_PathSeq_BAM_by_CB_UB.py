import pysam
from collections import defaultdict


# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile(pathseq.input[0], mode="rb")

iter = pathseq_bam.fetch()
output = []
UMI_dict = defaultdict(list)
# seg is an AlignedSegment object
for seg in iter:
    # not all records will have the CB tag and the UB tag
    if seg.has_tag("CB") and seg.has_tag("UB"):
        if (seg.get_tag(tag="CB") == snakemake.wildcards["cell"]):
            UMI_dict[seg.get_tag(tag="UB")].append(seg)


barcode_bam = pysam.AlignmentFile(snakemake.output[0], mode="wb", template=pathseq_bam)
for UMI in UMI_dict:
    # keep one read per UMI - the read with the highest mapping quality
    UMI_reads = UMI_dict[UMI]
    UMI_read = UMI_reads[0]
    for read in UMI_reads:
        if read.mapping_quality > UMI_read.mapping_quality:
            UMI_read = read
    barcode_bam.write(UMI_read)
barcode_bam.close()

pathseq_bam.close()
