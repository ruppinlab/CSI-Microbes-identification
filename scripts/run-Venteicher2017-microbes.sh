cd Venteicher2017/identify-microbes-workflow
sbatch \
--time=7-00:00:00 \
--cpus-per-task=8 \
--mem=8g \
--partition=norm,ccr \
scripts/run-snakemake.sh
