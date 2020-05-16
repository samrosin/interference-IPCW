#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=750m
#SBATCH -n 1
#SBATCH -t 00:10:00
module load r 

R CMD BATCH --no-save --no-restore dataset_gen_ni50_m500.R dataset_gen_ni50_m500-$SLURM_ARRAY_TASK_ID.Rout