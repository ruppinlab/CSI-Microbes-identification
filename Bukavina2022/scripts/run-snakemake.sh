#!/bin/bash

module load snakemake
snakemake \
--use-conda \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--group-components split_PathSeq_BAM_by_CB_UB=500 PathSeqScoreSpark=200 \
--local-cores 4 all
