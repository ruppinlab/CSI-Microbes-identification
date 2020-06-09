cd Chung2017
sbatch \
--time=7-00:00:00 \
--cpus-per-task=2 \
--mem=4g \
--partition=norm,ccr \
scripts/run-many-snakejobs.sh
