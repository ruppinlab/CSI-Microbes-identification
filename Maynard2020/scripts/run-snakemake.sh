#!/bin/bash

module load snakemake/6.0.5
snakemake \
--use-conda \
--nolock \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 100 \
--latency-wait 60 \
--keep-going \
--group-components score_PathSeq_cell_BAM=50 extract_unpaired_reads=100 sort_paired_reads=100 extract_paired_reads=100 split_PathSeq_BAM_by_RG=50 \
--local-cores 2 all
