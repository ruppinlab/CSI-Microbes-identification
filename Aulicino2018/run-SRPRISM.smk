# include Snakefile
include: "Snakefile"

# URLs
LT2_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"
D23580_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/538/085/GCF_900538085.1_D23580_liv/GCF_900538085.1_D23580_liv_genomic.fna.gz"
LT2_GFF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz"
D23580_GFF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/538/085/GCF_900538085.1_D23580_liv/GCF_900538085.1_D23580_liv_genomic.gff.gz"

wildcard_constraints:
    genome="LT2|D23580",
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

# SRPRISM_FILTER_COUNT_FILE = join("output", "SRPRISM", "{patient}", "{sample}", "{filter}_read_counts.tsv")
# SRPRISM_FILTER_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.{filter}.primary.sorted.bam")
# SRPRISM_NON_FILTER_BAM = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.non.{filter}.primary.sorted.bam")
# MPILEUP_FILE = join("output", "SRPRISM", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired.primary.sorted.mpileup")

GFF_READCOUNT_FILE = join("output", "{patient}", "{sample}-{plate}-{cell}", "{genome}-paired-count.gff")



localrules: intersect_BAM_GFF

rule all:
    input:
        expand(GFF_READCOUNT_FILE, zip, patient=cells["patient"], sample=cells["sample"], plate=cells["plate"], cell=cells["cell"], genome=["LT2"]*cells.shape[0]),
        expand(GFF_READCOUNT_FILE, zip, patient=cells["patient"], sample=cells["sample"], plate=cells["plate"], cell=cells["cell"], genome=["D23580"]*cells.shape[0]),
        expand(SRPRISM_COUNT_FILE, patient="Pt0", sample="S0", genome=["LT2", "D23580"]),

### rules for analyzing SRPRISM output ###
# filter only reads from cell of interest from the STAR bam file
rule split_STAR_unaligned_BAM_by_RG:
    group:
        "SRPRISM"
    input:
        join("output", "star", "{patient}-{sample}-{plate}", "_STARpe", "unaligned.bam"),
    output:
        join("output", "star", "{patient}-{sample}-{plate}-{cell}", "_STARpe", "unaligned.bam"),
    shell:
        "module load samtools && "
        "samtools view -h -b -r {wildcards.cell} {input} > {output}"

rule extract_FQ_files_from_BAM:
    group:
        "SRPRISM"
    input:
        join("output", "star", "{patient}-{sample}-{plate}-{cell}", "_STARpe", "unaligned.bam"),
    output:
        SRPRISM_INPUT_FQ1,
        SRPRISM_INPUT_FQ2
    shell:
        "module load bedtools && "
        "bamToFastq -i {input} -fq {output[0]} -fq2 {output[1]}"

rule intersect_BAM_GFF:
    input:
        SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM,
        GENOME_GFF
    output:
        GFF_READCOUNT_FILE
    shell:
        "module load bedtools && "
        "bedtools intersect -a {input[1]} -b {input[0]} -c > {output}"


# rule get_genome_coverage:
#     input:
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM,
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI,
#     shell:
#         "module load bedtools && "
#         "bedtools genomecov -d -ibam {input[0]} "



# rules for manipulating SPRISM bam files
# rule filter_reads:
#     input:
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAM,
#         SRPRISM_PROPER_PAIRED_PRIMARY_SORTED_BAI,
#         BED_FILTER_FILE
#     output:
#         SRPRISM_FILTER_BAM,
#         SRPRISM_NON_FILTER_BAM
#     shell:
#         "module load samtools && "
#         "samtools view -h -b -L {input[2]} -U {output[1]} {input[0]} > {output[0]}"

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


# def get_filtered_bam_files(wildcards):
#     return expand(SRPRISM_FILTER_BAM,  zip, patient=cells["patient"],
#            sample=cells["sample"], plate=cells["plate"], cell=cells["cell"],
#            genome=cells["infection"], filter=[wildcards.filter]*cells.shape[0])

# rule count_filtered_reads:
#     params:
#         exposed_cells["cell"]
#     conda:
#         "../envs/pysam-env.yaml"
#     input:
#         get_filtered_bam_files
#     output:
#         SRPRISM_FILTER_COUNT_FILE
#     script:
#         "../src/count_nreads.py"


# rule extract_filter_BED:
#     wildcard_constraints:
#         filter="rRNA|16S|protein_coding"
#     conda:
#         "../envs/pysam-env.yaml"
#     input:
#         GENOME_GTF
#     output:
#         BED_FILTER_FILE,
#     script:
#         "src/extract_rRNA_coordinates.py"



### Rules and functions for downloading Salmonella genome and transcriptome files ###
rule download_D23580_genome:
    wildcard_constraints:
        genome="D23580"
    params:
        D23580_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_LT2_genome:
    wildcard_constraints:
        genome="LT2"
    params:
        LT2_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_D23580_GFF:
    wildcard_constraints:
        genome="D23580"
    params:
        D23580_GFF_URL
    output:
        GENOME_GFF
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_LT2_GFF:
    wildcard_constraints:
        genome="LT2"
    params:
        LT2_GFF_URL
    output:
        GENOME_GFF
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"
