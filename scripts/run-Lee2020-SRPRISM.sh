cd Lee2020
sbatch \
--time=8:00:00 \
--cpus-per-task=2 \
--mem=4g \
--partition=norm,ccr \
scripts/run-SRPRISM.sh
