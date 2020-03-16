cd Yost2019/identify-microbes-workflow
sbatch \
--time=7-00:00:00 \
--cpus-per-task=16 \
--mem=8g \
--partition=norm,ccr \
scripts/run-snakemake.sh
