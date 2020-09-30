cd Puram2017
sbatch \
--time=7-00:00:00 \
--cpus-per-task=16 \
--mem=32g \
--partition=norm,ccr \
scripts/run-snakemake.sh
