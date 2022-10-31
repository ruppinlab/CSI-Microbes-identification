#!/bin/bash

module load snakemake
snakemake \
--use-conda \
--nolock \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 10 \
--latency-wait 60 \
--keep-going \
--group-components score_PathSeq_cell_BAM=50 extract_unpaired_reads=100 sort_paired_reads=100 extract_paired_reads=100 split_PathSeq_BAM_by_RG=400 split_STAR_unaligned_BAM_by_RG=50 extract_FQ_files_from_BAM=400 intersect_BAM_GFF=400 index_bam=400 sort_bam=400 convert_to_bam=400 exclude_non_proper_pairs=400 extract_primary_alignment=400 map_SRPRISM_genome_paired=400 \
--local-cores 4 \
--snakefile run-SRPRISM.smk
