#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=4g
#SBATCH -n 1
#SBATCH -t 08:00:00
module load r 

R CMD BATCH --no-save --no-restore results_gen_ni50_m500.R results_gen_ni50_m500-$SLURM_ARRAY_TASK_ID.Rout