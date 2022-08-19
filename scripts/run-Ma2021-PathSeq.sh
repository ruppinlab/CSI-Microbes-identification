cd Ma2021
sbatch \
--time=2-00:00:00 \
--cpus-per-task=4 \
--mem=8g \
--partition=norm,ccr \
scripts/run-snakemake.sh
