include: "Snakefile"
include: "../pathogen-discovery-rules/rules/PathSeq-10x.smk"


# PathSeq files
PATHSEQ_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam")
PATHSEQ_CELL_SCORE = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq.txt")


rule all:
    input:
        expand("output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.genus.csv", zip, patient=samples["patient"], sample=samples["sample"])


rule run_invade_seq:
    conda:
        "../envs/pysam-env.yaml"
    input:
        CR_bam=CR_BAM_FILE,
        barcodes=CR_RAW_BARCODES_FILE,
        pathseq_bam=PATHSEQ_FILTERED_BAM,
        pathseq_scores=join("output", "PathSeq", "{patient}-{sample}", "pathseq.txt"),
    output:
        read_name_pathseq="output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.readname",
        unmap_cbub_bam_file = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.bam",
        unmap_cbub_fasta_file = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.fasta",
        out_cell_list = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.unmap_cbub.list",
        out_readname_cell_path = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.raw.filtered_matrix.readnamepath",
        out_genus_file = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.genus.cell",
        output_UMI_table_csv = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.genus.csv",
        output_UMI_validate_table_csv  = "output/INVADEseq_CSI-Microbes-comparison/{patient}-{sample}.gex.filtered_matrix.validate.cell"
    shell:
        "python ../../Galeano-Nino-Bullman-Intratumoral-Microbiota_2022/patient_samples/INVADEseq.py "
        "{input.CR_bam} {wildcards.patient}-{wildcards.sample} {input.barcodes} {input.pathseq_bam} {input.pathseq_scores} "
        "{output.read_name_pathseq} {output.unmap_cbub_bam_file} {output.unmap_cbub_fasta_file} {output.out_cell_list} {output.out_readname_cell_path} {output.out_genus_file} {output.output_UMI_table_csv} {output.output_UMI_validate_table_csv}"
