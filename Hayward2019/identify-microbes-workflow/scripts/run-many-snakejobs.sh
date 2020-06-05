#!/bin/bash
INPUT=scripts/data.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read start stop
do
  sbatch \
  --time=7-00:00:00 \
  --cpus-per-task=2 \
  --mem=4g \
  --partition=norm,ccr \
  scripts/run-snakemake.sh $start $stop

done < $INPUT