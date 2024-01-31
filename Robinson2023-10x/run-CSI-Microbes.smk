include: "Snakefile"
include: "../pathogen-discovery-rules/rules/PathSeq-10x.smk"


# PathSeq files
PATHSEQ_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam")
PATHSEQ_CELL_SCORE = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq.txt")


rule all:
    input:
        expand(PATHSEQ_CELL_SCORE, zip, patient=cells["patient"], sample=cells["sample"], cell=cells["barcode"])

