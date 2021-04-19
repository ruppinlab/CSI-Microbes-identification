# include Snakefile
include: "Snakefile"

# URLs
SL1344_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"
SL1344_GTF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gtf.gz"

samples["genome"] = "SL1344"

include: "../RNA-snakemake-rules/rules/SRPRISM-unpaired.smk"

# SRPRISM files
GENOME_FA = join("raw", "{genome}.fa")
GENOME_GTF = join("raw", "{genome}.gtf")
TRIMMED_FQ1 = join("FASTQ", "unmapped", "trimmed", "{patient}-{sample}_1.fastq.gz")
SRPRISM_UNPAIRED_INPUT_FQ = join("output", "SRPRISM", "{patient}", "{sample}", "unaligned_3.fq")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.bam")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.bam.bai")
#SRPRISM_TAG_BAM = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.withtags.bam")
SRPRISM_READ_COUNT = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-nreads.tsv") # reads for a sample
SRPRISM_CB_UMI_TABLE = join("output", "SRPRISM", "{patient}", "{sample}", "all-CB-UMI-table-{genome}.tsv") # table of all mapped reads with CB and UMI for a sample
SRPRISM_CB_UMI_COUNT = join("output", "SRPRISM", "{patient}", "{sample}", "CB-UMI-count-{genome}.tsv") # table of the number of UMIs for each CB for a sample


rule SRPRISM_output:
    input:
        expand(SRPRISM_READ_COUNT, zip, patient=samples["patient"], sample=samples["sample"], genome=samples["genome"]),
        expand(SRPRISM_CB_UMI_TABLE, zip, patient=samples["patient"], sample=samples["sample"], genome=samples["genome"]),
        expand(SRPRISM_CB_UMI_COUNT, zip, patient=samples["patient"], sample=samples["sample"], genome=samples["genome"]),

rule add_CR_tags_SRPRISM:
    conda:
        "../envs/pysam-env.yaml"
    input:
        CR_BAM_FILE,
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM
    output:
        SRPRISM_CB_UMI_TABLE,
        SRPRISM_CB_UMI_COUNT
    script:
        "src/add_CR_tags_to_SRPRISM_bam.py"

# add the CB and UB tag from the CellRanger BAM to the SRPRISM BAM
# rule add_CB_UB_tags_to_SRPRISM_BAM:
#     conda:
#         "../pathogen-discovery-rules/envs/pysam.yaml"
#     input:
#         CR_BAM_FILE,
#         SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM,
#     output:
#         SRPRISM_TAG_BAM,
#     script:
#         "../pathogen-discovery-rules/src/add_tags_to_PathSeq_bam.py"

# rule to count the total number of reads from a BAM sample
rule count_total_n_reads:
    input:
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM
    output:
        SRPRISM_READ_COUNT
    shell:
        "module load samtools && "
        "samtools view -c {input} > {output}"

rule move_FQ_for_SRPRISM:
    input:
        TRIMMED_FQ1
    output:
        temp(SRPRISM_UNPAIRED_INPUT_FQ)
    shell:
        "zcat -c {input} > {output}"



### rules to download SL1344 reference files ###

rule download_SL1344_genome:
    wildcard_constraints:
        genome="SL1344"
    params:
        SL1344_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_SL1344_GTF:
    wildcard_constraints:
        genome="SL1344"
    params:
        SL1344_GTF_URL
    output:
        GENOME_GTF
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"
