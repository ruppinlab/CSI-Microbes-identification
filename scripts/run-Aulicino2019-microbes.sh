cd Aulicino2018/identify-microbes-workflow
sbatch \
--time=7-00:00:00 \
--mem=4g \
--partition=norm,ccr \
scripts/run-snakemake.sh
