#!/bin/bash
INPUT=scripts/data.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read batch
do
  sbatch \
  --time=7-00:00:00 \
  --cpus-per-task=16 \
  --mem=16g \
  --partition=norm,ccr \
  scripts/run-snakemake.sh $batch

done < $INPUT
