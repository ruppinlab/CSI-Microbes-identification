cd build-PathSeq-microbe-files
sbatch \
--time=4-00:00:00 \
--cpus-per-task=4 \
--mem=4g \
--partition=norm,ccr \
scripts/run-snakemake.sh
