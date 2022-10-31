# include Snakefile
include: "Snakefile"

# URLs
Fn_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/457/555/GCF_001457555.1_NCTC10562/GCF_001457555.1_NCTC10562_genomic.fna.gz"
Fn_GFF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/457/555/GCF_001457555.1_NCTC10562/GCF_001457555.1_NCTC10562_genomic.gff.gz"

wildcard_constraints:
    genome="Fn",
    filter="rRNA|16S|protein_coding"

# include rules
include: "../RNA-snakemake-rules/rules/SRPRISM-paired.smk"

# Intermediate Files for SRPRISM
GENOME_FA = join("raw", "{genome}.fa")
GENOME_FAI = join("raw", "{genome}.fa.fai")
GENOME_GFF = join("raw", "{genome}.gff")
BED_FILTER_FILE=join("raw", "{genome}_{filter}.bed")

# SRPRISM Output Files
SRPRISM_INPUT_FQ1 = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "unaligned_1.fq")
SRPRISM_INPUT_FQ2 = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "unaligned_2.fq")
SRPRISM_PROPER_PAIRED_PRIMARY_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-proper-paired.primary.bam")
SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.primary.sorted.bam")
SRPRISM_COUNT_FILE = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}_read_counts.tsv")

GFF_READCOUNT_FILE = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")




rule SRPRISM_files:
    input:
        expand(GFF_READCOUNT_FILE, zip, patient=cells["patient"], sample=cells["sample"], plate=cells["plate"], cell=cells["cell"], genome=["Fn"]*cells.shape[0]),
        expand(SRPRISM_COUNT_FILE, patient="Pt0", sample="S0", genome=["Fn"]),

### rules for analyzing SRPRISM output ###
# filter only reads from cell of interest from the STAR bam file
rule split_STAR_unaligned_BAM_by_RG:
    group:
        "split_STAR_unaligned_BAM_by_RG"
    input:
        join("output", "star", "{patient}-{sample}-{plate}", "_STARpe", "unaligned.bam"),
    output:
        join("output", "star", "{patient}-{sample}-{plate}-{cell}", "_STARpe", "unaligned.bam"),
    shell:
        "module load samtools && "
        "samtools view -h -b -r {wildcards.cell} {input} > {output}"

rule extract_FQ_files_from_BAM:
    group:
        "extract_FQ_files_from_BAM"
    input:
        join("output", "star", "{patient}-{sample}-{plate}-{cell}", "_STARpe", "unaligned.bam"),
    output:
        SRPRISM_INPUT_FQ1,
        SRPRISM_INPUT_FQ2
    shell:
        "module load bedtools && "
        "bamToFastq -i {input} -fq {output[0]} -fq2 {output[1]}"

rule intersect_BAM_GFF:
    group:
        "intersect_BAM_GFF"
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM,
        GENOME_GFF
    output:
        GFF_READCOUNT_FILE
    shell:
        "module load bedtools && "
        "bedtools intersect -a {input[1]} -b {input[0]} -c > {output}"

# count the reads and combine
rule count_total_reads:
    params:
        cells["cell"]
    conda:
        "../envs/pysam-env.yaml"
    input:
        expand(SRPRISM_PROPER_PAIRED_PRIMARY_BAM,  zip, patient=cells["patient"],
               sample=cells["sample"], plate=cells["plate"], cell=cells["cell"],
               genome=["{genome}"]*cells.shape[0])
    output:
        SRPRISM_COUNT_FILE
    script:
        "../src/count_nreads.py"

### rules to download Fn reference files ###
rule download_Fn_genome:
    wildcard_constraints:
        genome="Fn"
    params:
        Fn_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_Fn_GFF:
    wildcard_constraints:
        genome="Fn"
    params:
        Fn_GFF_URL
    output:
        GENOME_GFF
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"
