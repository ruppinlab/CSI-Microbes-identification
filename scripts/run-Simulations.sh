cd Simulations
sbatch \
--time=4-00:00:00 \
--cpus-per-task=2 \
--mem=4g \
--partition=norm,ccr \
scripts/run-snakemake.sh
