cd Lee2020
sbatch \
--time=7-00:00:00 \
--cpus-per-task=16 \
--mem=48g \
--partition=norm,ccr \
scripts/run-snakemake.sh
