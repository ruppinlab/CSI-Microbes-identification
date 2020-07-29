#!/bin/bash
INPUT=scripts/data.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read patient plate
do
  sbatch \
  --time=7-00:00:00 \
  --cpus-per-task=8 \
  --mem=16g \
  --partition=norm,ccr \
  scripts/run-snakemake.sh $patient $plate

done < $INPUT
