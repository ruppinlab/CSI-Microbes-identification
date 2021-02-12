cd test-10x 
sbatch \
--time=2-00:00:00 \
--cpus-per-task=16 \
--mem=64g \
--partition=norm,ccr \
scripts/run-snakemake.sh
