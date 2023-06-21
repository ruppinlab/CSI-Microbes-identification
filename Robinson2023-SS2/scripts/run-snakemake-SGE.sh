#!/bin/bash --login

conda activate snakemake

snakemake \
--profile sge \
--local-cores 4 all \
&> stdout_snakemake.log
