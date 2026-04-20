#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --chdir=/data/WIPH-CanDetect/HealthEco/scripts/job_scripts/
#SBATCH --output=/data/WIPH-CanDetect/HealthEco/scripts/job_scripts/run_records/%A_%a.out
#SBATCH --array=1-16:1  

iter=${SLURM_ARRAY_TASK_ID}
echo ${iter}

module load R/4.5.1

# Replace the following line with a program or command
Rscript master_tp_optimise_by2Char_4hpc.R ${iter} "/data/WIPH-CanDetect/HealthEco/output/optim/" "age70plus" "female" "RS" "500000" "opt_upd2026" 