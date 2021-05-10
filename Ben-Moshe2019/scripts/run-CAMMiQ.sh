#!/bin/bash

module load snakemake/6.0.5
snakemake \
--use-conda \
--nolock \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 10 \
--latency-wait 60 \
--keep-going \
--group-components split_PathSeq_BAM_by_CB_UB=500 PathSeqScoreSpark=500 convert_BAMs_to_fastq_dir=500 \
--local-cores 4 \
--snakefile run-CAMMiQ.smk
