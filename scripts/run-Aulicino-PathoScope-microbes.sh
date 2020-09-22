cd Aulicino-PathoScope
sbatch \
--time=1-00:00:00 \
--cpus-per-task=2 \
--mem=4g \
--partition=norm,ccr \
scripts/run-snakemake.sh
