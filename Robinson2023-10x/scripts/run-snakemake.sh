#!/bin/bash --login

conda activate snakemake
ml Singularity/3.6.4
snakemake \
--use-conda \
--use-singularity \
--rerun-incomplete \
--cluster-config config/crick-cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} " \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--group-components split_PathSeq_BAM_by_CB_UB=10000 PathSeqScoreSpark=10000 \
--local-cores 4 all
