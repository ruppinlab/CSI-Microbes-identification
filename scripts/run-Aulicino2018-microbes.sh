cd Aulicino2018/identify-microbes-workflow
sbatch \
--time=3-00:00:00 \
--cpus-per-task=16 \
--mem=8g \
--partition=norm,ccr \
scripts/run-snakemake.sh
