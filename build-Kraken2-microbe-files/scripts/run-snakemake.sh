#!/bin/bash --login

conda activate snakemake
ml Singularity/3.6.4
snakemake \
--use-conda \
--use-singularity \
--singularity-args "--bind data:/data" \
--rerun-incomplete \
--cluster-config config/crick-cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--local-cores 4 build_kraken_db 
