cd Simmons-Collaboration
sbatch \
--time=7-00:00:00 \
--cpus-per-task=24 \
--mem=72g \
--partition=norm,ccr \
scripts/run-snakemake.sh
