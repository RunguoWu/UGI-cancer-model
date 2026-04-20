#!/bin/bash
#$ -pe smp 8               
#$ -l h_vmem=4G       
#$ -l h_rt=24:0:0 
        
#$ -cwd                       
#$ -j y
#$ -o /data/WIPH-CanDetect/HealthEco/scripts/index_1yr/run_records/                      
#$ -t 1-16:1     

iter=${SGE_TASK_ID}

echo ${iter}

module load R/4.5.1

# Replace the following line with a program or command
Rscript master_tp_optimise_by2Char_4hpc_1yrIndex.R ${iter} "/data/WIPH-CanDetect/HealthEco/output/optim/" "age70plus" "female" "RS" "500000" "opt_1yrIndex" 