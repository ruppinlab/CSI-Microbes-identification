cd Lambrechts2018
sbatch \
--time=2-00:00:00 \
--cpus-per-task=2 \
--mem=16g \
--partition=norm,ccr \
scripts/run-snakemake.sh
