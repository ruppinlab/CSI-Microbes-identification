sbatch \
--time=7-00:00:00 \
--cpus-per-task=4 \
--mem=16g \
--partition=norm,ccr \
scripts/run-snakemake.sh
