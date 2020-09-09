cd Ben-Moshe2019
sbatch \
--time=7-00:00:00 \
--cpus-per-task=64 \
--mem=96g \
--partition=norm,ccr \
scripts/run-snakemake.sh
