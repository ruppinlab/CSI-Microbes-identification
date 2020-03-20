#!/bin/bash

module load snakemake
snakemake \
--use-conda \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition=norm,ccr --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--local-cores 16 all
