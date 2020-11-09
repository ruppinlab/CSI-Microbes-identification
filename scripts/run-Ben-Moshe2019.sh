cd Ben-Moshe2019
sbatch \
--time=1-00:00:00 \
--cpus-per-task=4 \
--mem=10g \
--partition=norm,ccr \
scripts/run-snakemake.sh
