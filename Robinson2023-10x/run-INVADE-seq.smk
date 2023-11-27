include: "Snakefile"


rule all:
    expand("output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.genus.csv", patient="P1", sample=["SCAF2961_1_Uninfected", "SCAF2962_2_HK", "SCAF2963_3_Live", "SCAF2965_5_Live"]),

### rules to run PathSeq and get per-cell microbial read abundances ###
rule PathSeqPipelineSpark:
    input:
        bam_file = CR_BAM_FILE,
        host_bwa_image = config["PathSeq"]["host_img"],
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_dict_file = config["PathSeq"]["microbe_dict"],
        host_hss_file = config["PathSeq"]["host_bfi"],
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    output:
        pathseq_bam = join("output", "INVADEseq" "PathSeq", "{patient}-{sample}", "pathseq.bam"),
        pathseq_output = join("output", "INVADEseq", "PathSeq", "{patient}-{sample}", "pathseq.txt"),
        filter_metrics = join("output", "INVADEseq", "PathSeq", "{patient}-{sample}", "filter-metrics.txt"),
        score_metrics = join("output", "INVADEseq", "PathSeq", "{patient}-{sample}", "score-metrics.txt"),
    benchmark:
        "benchmarks/{patient}-{sample}.PathSeqPipelineSpark_host_filter_single.txt"
    run:
    shell(
            "ml GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8 && "
            "gatk PathSeqPipelineSpark "
            "--filter-duplicates false "
            "--min-score-identity .7 "
            "--input '{input.bam_file}' "
            "--is-host-aligned false "
            "--min-clipped-read-length 60 "
            "--filter-bwa-image /lscratch/$SLURM_JOBID/{input.host_bwa_image} "
            "--kmer-file /lscratch/$SLURM_JOBID/{input.host_hss_file} "
            "--microbe-bwa-image /lscratch/$SLURM_JOBID/{input.microbe_bwa_image} "
            "--microbe-dict /lscratch/$SLURM_JOBID/{input.microbe_dict_file} "
            "--taxonomy-file /lscratch/$SLURM_JOBID/{input.taxonomy_db} "
            "--output '{output.pathseq_bam}' "
            "--scores-output '{output.pathseq_output}' "
            "--filter-metrics '{output.filter_metrics}' "
            "--score-metrics '{output.score_metrics}' "
            '--java-options "-Xmx96g -Xms96G -XX:+UseG1GC -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" '
            '--spark-master local[8] '
            + config["params"]["PathSeq"]
        )



rule run_invade_seq:
    conda:
        "../envs/pysam-env.yaml"
    input:
        CR_bam=CR_BAM_FILE,
        barcodes=CR_BARCODES_FILE,
        pathseq_bam=join("output", "INVADEseq" "PathSeq", "{patient}-{sample}", "pathseq.bam"),
        pathseq_scores=join("output", "INVADEseq", "PathSeq", "{patient}-{sample}", "pathseq.txt"),
    output:
        read_name_pathseq="output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.readname",
        unmap_cbub_bam_file = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.bam",
        unmap_cbub_fasta_file = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.fasta",
        out_cell_list = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.list",
        out_readname_cell_path = "output/INVADEseq_raw/{patient}-{sample}.gex.raw.filtered_matrix.readnamepath",
        out_genus_file = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.genus.cell",
        output_UMI_table_csv = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.genus.csv",
        output_UMI_validate_table_csv  = "output/INVADEseq_raw/{patient}-{sample}.gex.filtered_matrix.validate.cell"
    shell:
        "python ../../Galeano-Nino-Bullman-Intratumoral-Microbiota_2022/patient_samples/INVADEseq.py "
        "{input.CR_bam} {wildcards.patient}-{wildcards.sample} {input.barcodes} {input.pathseq_bam} {input.pathseq_scores} "
        "{output.read_name_pathseq} {output.unmap_cbub_bam_file} {output.unmap_cbub_fasta_file} {output.out_cell_list} {output.out_readname_cell_path} {output.out_genus_file} {output.output_UMI_table_csv} {output.output_UMI_validate_table_csv}"
