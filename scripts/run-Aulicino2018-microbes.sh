cd Aulicino2018
sbatch \
--time=2-00:00:00 \
--cpus-per-task=2 \
--mem=4g \
--partition=norm,ccr \
scripts/run-many-snakejobs.sh
