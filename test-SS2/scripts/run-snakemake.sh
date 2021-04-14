#!/bin/bash

module load snakemake/6.0.5
snakemake \
--use-conda \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--group-components download_FASTQ_from_SRA=15 compress_FASTQ_File=15 run_fastp=15 split_PathSeq_BAM_by_RG=15 extract_paired_reads=15 sort_paired_reads=15 extract_unpaired_reads=15 score_PathSeq_cell_BAM=15 \
--local-cores 4 all
