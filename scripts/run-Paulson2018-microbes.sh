cd Paulson2018
sbatch \
--time=7-00:00:00 \
--cpus-per-task=1 \
--mem=4g \
--partition=norm,ccr \
scripts/run-snakemake.sh
