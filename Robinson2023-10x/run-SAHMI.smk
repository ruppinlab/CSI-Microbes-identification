include: "Snakefile"


localrules: clean_kraken_output, concatenate_fq_gz_files


rule all:
    input:
        expand("output/SAHMI/{patient}/{sample}.counts.txt", patient="P1", sample=["SCAF2961_1_Uninfected", "SCAF2962_2_HK", "SCAF2963_3_Live", "SCAF2965_5_Live"])


rule concatenate_fq_gz_files:
    input:
        unpack(get_sra_fq_files_by_sample)
    output:
        CellBarcode_fq = CB_SAMPLE_FASTQ_FILE,
        cDNA_fq = cDNA_SAMPLE_FASTQ_FILE, 
    shell:
        "cat {input.CellBarcode_fq} > {output.CellBarcode_fq} && "
        "cat {input.cDNA_fq} > {output.cDNA_fq}"


rule run_kraken:
    singularity:
        "library://wir963/csi-microbes/kraken2" 
    params:
        "output/kraken2/{patient}/{sample}#.fq"
    input:
        CB_SAMPLE_FASTQ_FILE,
        cDNA_SAMPLE_FASTQ_FILE 
    output:
        fq1="output/kraken2/{patient}/{sample}_1.fq",
        fq2="output/kraken2/{patient}/{sample}_2.fq",
        output="output/kraken2/{patient}/{sample}.kraken.output.txt",
        report="output/kraken2/{patient}/{sample}.kraken.report.txt"
    shell:
        "export OMP_NUM_THREADS=16 && "
        "kraken2 --db microbev1 "
        "--paired "
        "--threads 16 "
        "--use-names "
        "--report-minimizer-data "
        "--classified-out {params} "
        "--output {output.output} "
        "--report {output.report} "
        "{input} "


rule clean_kraken_output:
    input:
        "output/kraken2/{patient}/{sample}.kraken.report.txt"
    output:
        std="output/kraken2/{patient}/{sample}.kraken.report.std.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt"
    shell:
        "cut -f1-3,6-8 {input} > {output.std} && "
        "../../SAHMI/functions/kreport2mpa.py -r {output.std} -o {output.mpa} --intermediate-ranks"


rule extract_microbiome_reads_fq1:
    conda:
        "../envs/SAHMI-env.yaml"
    params:
        "output/SAHMI/{patient}/r1/"
    input:
        fq="output/kraken2/{patient}/{sample}_1.fq",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    output:
        "output/SAHMI/{patient}/r1/{sample}.fa"
    shell:
        "Rscript ../../SAHMI/functions/extract_microbiome_reads.r --sample_name {wildcards.sample} "
        "--fq {input.fq} --kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params}"


rule extract_microbiome_reads_fq2:
    conda:
        "../envs/SAHMI-env.yaml"
    params:
        "output/SAHMI/{patient}/r2/"
    input:
        fq="output/kraken2/{patient}/{sample}_2.fq",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    output:
        "output/SAHMI/{patient}/r2/{sample}.fa"
    shell:
        "Rscript ../../SAHMI/functions/extract_microbiome_reads.r --sample_name {wildcards.sample} "
        "--fq {input.fq} --kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params}"


rule extract_microbiome_output:
    conda:
        "../envs/SAHMI-env.yaml"
    params:
        "output/SAHMI/{patient}/"
    input:
        out="output/kraken2/{patient}/{sample}.kraken.output.txt",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    output:
        "output/SAHMI/{patient}/{sample}.microbiome.output.txt"
    shell:
        "Rscript ../../SAHMI/functions/extract_microbiome_output.r --sample_name {wildcards.sample} "
        "--output_file {input.out} --kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params}"

# need to run sckmer paired and fa1 must be the read containing the cell barcodes; annoying this isn't specified in the documentation
rule sckmer_paired_3p_v3:
    conda:
        "../envs/SAHMI-env.yaml"
    wildcard_constraints:
        sample="SCAF2961_1_Uninfected|SCAF2962_2_HK|SCAF2963_3_Live"
    params:
        "output/SAHMI/{patient}/"
    input:
        fa1="output/SAHMI/{patient}/r1/{sample}.fa",
        fa2="output/SAHMI/{patient}/r2/{sample}.fa",
        micro="output/SAHMI/{patient}/{sample}.microbiome.output.txt",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    output:
        "output/SAHMI/{patient}/{sample}.sckmer.txt"
    shell:
        "Rscript ../../SAHMI/functions/sckmer.r --sample_name {wildcards.sample} --fa1 {input.fa1} --fa2 {input.fa2} --microbiome_output_file {input.micro} "
        "--kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params} --cb_len 16 --umi_len 12 "

rule sckmer_paired_5p:
    conda:
        "../envs/SAHMI-env.yaml"
    wildcard_constraints:
        sample="SCAF2965_5_Live"
    params:
        "output/SAHMI/{patient}/"
    input:
        fa1="output/SAHMI/{patient}/r1/{sample}.fa",
        fa2="output/SAHMI/{patient}/r2/{sample}.fa",
        micro="output/SAHMI/{patient}/{sample}.microbiome.output.txt",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    output:
        "output/SAHMI/{patient}/{sample}.sckmer.txt"
    shell:
        "Rscript ../../SAHMI/functions/sckmer.r --sample_name {wildcards.sample} --fa1 {input.fa1} --fa2 {input.fa2} --microbiome_output_file {input.micro} "
        "--kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params} --cb_len 16 --umi_len 10 "

rule barcode_denoising:
    conda:
        "../envs/SAHMI-env.yaml"
    input:
        kmer="output/SAHMI/{patient}/{sample}.sckmer.txt",
        report="output/kraken2/{patient}/{sample}.kraken.report.txt"
    output:
        hits="output/SAHMI/{patient}/{sample}.barcode.kmer.hits.tsv",
        rpmm="output/SAHMI/{patient}/{sample}.kraken.report.rpmm.tsv"
    script:
        "src/barcode_denoising.R" 

rule cell_line_quantile_test:
    input:
        kmer="output/SAHMI/{patient}/{sample}.barcode.kmer.hits.tsv",
        rpmm="output/SAHMI/{patient}/{sample}.kraken.report.rpmm.tsv",
        cell_line="../../SAHMI/Table S4.xlsx"
    output:
        hits="output/SAHMI/{patient}/{sample}.cell_line_quantile_hits.tsv",
        taxa="output/SAHMI/{patient}/{sample}.cell_line_quantile_hits_taxa.tsv"
    script:
        "src/cell_line_quantile_test.py"

rule taxa_counts_3p_v3:
    conda:
        "../envs/SAHMI-env.yaml"
    wildcard_constraints:
        sample="SCAF2961_1_Uninfected|SCAF2962_2_HK|SCAF2963_3_Live"
    input:
        fa1="output/SAHMI/{patient}/r1/{sample}.fa",
        fa2="output/SAHMI/{patient}/r2/{sample}.fa",
        taxa="output/SAHMI/{patient}/{sample}.cell_line_quantile_hits_taxa.tsv",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    params:
        "output/SAHMI/{patient}/"
    output:
        "output/SAHMI/{patient}/{sample}.counts.txt"
    shell:
        "Rscript ../../SAHMI/functions/taxa_counts.r --sample_name {wildcards.sample} --fa1 {input.fa1} --fa2 {input.fa2} "
        "--taxa {input.taxa} --kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params} --cb_len 16 --umi_len 12 "

rule taxa_counts_5p:
    conda:
        "../envs/SAHMI-env.yaml"
    wildcard_constraints:
        sample="SCAF2965_5_Live"
    input:
        fa1="output/SAHMI/{patient}/r1/{sample}.fa",
        fa2="output/SAHMI/{patient}/r2/{sample}.fa",
        taxa="output/SAHMI/{patient}/{sample}.cell_line_quantile_hits_taxa.tsv",
        rep="output/kraken2/{patient}/{sample}.kraken.report.txt",
        mpa="output/kraken2/{patient}/{sample}.kraken.report.mpa.txt",
    params:
        "output/SAHMI/{patient}/"
    output:
        "output/SAHMI/{patient}/{sample}.counts.txt"
    shell:
        "Rscript ../../SAHMI/functions/taxa_counts.r --sample_name {wildcards.sample} --fa1 {input.fa1} --fa2 {input.fa2} "
        "--taxa {input.taxa} --kraken_report {input.rep} --mpa_report {input.mpa} --out_path {params} --cb_len 16 --umi_len 10 "
