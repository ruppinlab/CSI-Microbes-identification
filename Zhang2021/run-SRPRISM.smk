# include
include: "Snakefile"

# URLs
Ps_genome_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/343/155/GCF_013343155.1_ASM1334315v1/GCF_013343155.1_ASM1334315v1_genomic.fna.gz"
Ps_GFF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/343/155/GCF_013343155.1_ASM1334315v1/GCF_013343155.1_ASM1334315v1_genomic.gff.gz"

include: "../RNA-snakemake-rules/rules/SRPRISM-unpaired.smk"

# SRPRISM files
GENOME_FA = join("raw", "{genome}.fa")
GENOME_GFF = join("raw", "{genome}.gff")
TRIMMED_FQ1 = join("FASTQ", "unmapped", "trimmed", "{patient}-{sample}_1.fastq.gz")
SRPRISM_UNPAIRED_INPUT_FQ = join("output", "SRPRISM", "{patient}", "{sample}", "unaligned_3.fq")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.bam")
SRPRISM_UNPAIRED_PRIMARY_SORTED_BAI = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.bam.bai")
#SRPRISM_TAG_BAM = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired.primary.sorted.withtags.bam")
SRPRISM_CB_UMI_TABLE = join("output", "SRPRISM", "{patient}", "{sample}", "all-CB-UMI-table-{genome}.tsv") # table of all mapped reads with CB and UMI for a sample
SRPRISM_CB_UMI_COUNT = join("output", "SRPRISM", "{patient}", "{sample}", "CB-UMI-count-{genome}.tsv") # table of the number of UMIs for each CB for a sample

GFF_READCOUNT_FILE = join("output", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired-count.gff")

localrules: download_Ps_genome, download_Ps_GFF, move_FQ_for_SRPRISM

samples = samples.loc[samples["patient"].isin(["P83", "P79", "P8"])]

rule SRPRISM_output:
    input:
        expand(GFF_READCOUNT_FILE, patient=samples["patient"], sample=samples["sample"], genome="Fn"),
        #expand(SRPRISM_CB_UMI_TABLE, zip, patient=samples["patient"], sample=samples["sample"], genome=samples["genome"]),
        #expand(SRPRISM_CB_UMI_COUNT, zip, patient=samples["patient"], sample=samples["sample"], genome=samples["genome"]),

# add cell barcode and UMI tags to the BAM and use to count UMIs per cell
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
        "../Ben-Moshe2019/src/add_CR_tags_to_SRPRISM_bam.py"

# get a read count per gene per sample file
rule intersect_BAM_GFF:
    group:
        "intersect_BAM_GFF"
    input:
        SRPRISM_UNPAIRED_PRIMARY_SORTED_BAM,
        GENOME_GFF
    output:
        GFF_READCOUNT_FILE
    shell:
        "module load bedtools && "
        "bedtools intersect -a {input[1]} -b {input[0]} -c > {output}"

rule move_FQ_for_SRPRISM:
    input:
        TRIMMED_FQ1
    output:
        temp(SRPRISM_UNPAIRED_INPUT_FQ)
    shell:
        "zcat -c {input} > {output}"





### rules to download Fn reference files ###
rule download_Ps_genome:
    wildcard_constraints:
        genome="Ps"
    params:
        Ps_genome_URL
    output:
        GENOME_FA
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"

rule download_Ps_GFF:
    wildcard_constraints:
        genome="Ps"
    params:
        Ps_GFF_URL
    output:
        GENOME_GFF
    shell:
        "wget -O - {params[0]} | gunzip -c > {output}"