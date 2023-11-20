#!/bin/bash --login

conda activate snakemake
snakemake \
--use-conda \
--nolock \
--rerun-incomplete \
--cluster-config config/crick-cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads}" \
--jobs 10 \
--latency-wait 60 \
--keep-going \
--group-components FASTQ=10 score_PathSeq_cell_BAM=25 extract_unpaired_reads=25 sort_paired_reads=25 extract_paired_reads=25 split_PathSeq_BAM_by_RG=25 \
--local-cores 4 all
