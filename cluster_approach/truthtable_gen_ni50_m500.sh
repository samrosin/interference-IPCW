#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=500m
#SBATCH -n 1
#SBATCH -t 00:20:00
module load r 

R CMD BATCH --no-save --no-restore truthtable_gen_ni50_m500.R truthtable_gen_ni50_m500.Rout