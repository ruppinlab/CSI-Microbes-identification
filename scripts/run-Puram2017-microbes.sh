cd Puram2017
sbatch \
--time=7-00:00:00 \
--cpus-per-task=32 \
--mem=48g \
--partition=norm,ccr \
scripts/run-snakemake.sh
