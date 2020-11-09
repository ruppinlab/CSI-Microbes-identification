import pysam
from collections import defaultdict
import pandas as pd

# load and index the CellRanger BAM file
cr_bam = pysam.AlignmentFile(snakemake.input[0], mode="rb")
cr_idx = pysam.IndexedReads(cr_bam)
cr_idx.build()

# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile(snakemake.input[1], mode="rb")

output = {} 
output["CB"] = []
output["UMI"] = []
# d = defaultdict(lambda: defaultdict(list))
# seg is an AlignedSegment object
for seg in pathseq_bam.fetch(until_eof=True):
    # returns an IteratorRowSelection object, which contains one or more AlignedSegment object
    cr_list = list(cr_idx.find(seg.query_name))
    # we assume that all records belonging to the same query name will have the same CB/UB tag
    # not all records will have the CB tag and the UB tag
    if cr_list[0].has_tag("CB") and cr_list[0].has_tag("UB"):
        output["CB"].append(cr_list[0].get_tag(tag="CB"))
        output["UMI"].append(cr_list[0].get_tag(tag="UB"))

print(output)
df = pd.DataFrame(data=output)
print(df)
df.to_csv(snakemake.output[0], sep="\t")
print(df.groupby("CB")["UMI"].nunique())
out_df = df.groupby("CB")["UMI"].nunique()
out_df.to_csv(snakemake.output[1], sep="\t")



cr_bam.close()
pathseq_bam.close()
