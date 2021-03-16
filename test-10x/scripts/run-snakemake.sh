#!/bin/bash

#module load snakemake
# temporarily use conda env for snakemake=6.0.5 because there is no module yet
# conda activate snakemake-env
snakemake \
--use-conda \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--group-components split_PathSeq_BAM_by_CB_UB=12 PathSeqScoreSpark=12 \
--local-cores 4 all
