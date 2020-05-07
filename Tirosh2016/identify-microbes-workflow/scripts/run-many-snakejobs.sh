#!/bin/bash
INPUT=scripts/test_data.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read patient, celltype, start, stop
do
  sbatch \
  --time=7-00:00:00 \
  --cpus-per-task=2 \
  --mem=4g \
  --partition=norm,ccr \
  scripts/run-snakemake.sh $patient $celltype $start $stop

done < $INPUT
