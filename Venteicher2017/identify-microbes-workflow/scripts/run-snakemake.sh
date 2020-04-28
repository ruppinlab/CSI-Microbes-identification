#!/bin/bash
INPUT=scripts/data.csv
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read patient plate
do
  module load snakemake
  snakemake \
  --use-conda \
  --nolock \
  --rerun-incomplete \
  --cluster-config config/cluster.json \
  --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
  --jobs 100 \
  --latency-wait 60 \
  --keep-going \
  --config patient=$patient plate=$plate \
  --local-cores 16 all
done < $INPUT
