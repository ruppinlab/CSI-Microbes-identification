cd Slyper2020
sbatch \
--time=7-00:00:00 \
--cpus-per-task=24 \
--mem=64g \
--partition=norm,ccr \
scripts/run-snakemake.sh
