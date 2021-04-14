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
--group-components download_FASTQ_from_SRA=15 compress_FASTQ_File=12 run_fastp=12 \
--local-cores 4 all
