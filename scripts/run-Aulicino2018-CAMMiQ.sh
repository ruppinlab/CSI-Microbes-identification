cd Aulicino2018
sbatch \
--time=1-00:00:00 \
--cpus-per-task=4 \
--mem=4g \
--partition=norm,ccr \
scripts/run-CAMMiQ.sh
