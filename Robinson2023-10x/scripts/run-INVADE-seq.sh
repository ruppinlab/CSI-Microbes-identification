#!/bin/bash --login

conda activate snakemake
snakemake \
--use-conda \
--rerun-incomplete \
--cluster-config config/biowulf-cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} " \
--jobs 100 \
--latency-wait 60 \
--keep-going \
-s run-INVADE-seq.smk \
--local-cores 4 all
