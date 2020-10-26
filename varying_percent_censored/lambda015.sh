#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=8g
#SBATCH -n 1
#SBATCH -t 26:30:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=srosin@live.unc.edu
module load r/4.0.1    

#define variables
RAN_SEED=2020
M=500
N_I=10
N_SIMS=1000
SUBDIR_NAME=lambda_015
LAMBDA_COMPONENT=.015

srun --output=/dev/null --error=/dev/null R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME $LAMBDA_COMPONENT" dataset_gen_lambda.R dataset_gen_$SUBDIR_NAME.Rout

sbatch --output=/dev/null --error=/dev/null --time=04:30:00 --mem=2g --array=1-$N_SIMS --job-name=results_gen_$SUBDIR_NAME --wait R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME $LAMBDA_COMPONENT" results_gen_lambda.R results_gen_$SUBDIR_NAME.Rout 

sbatch --output=/dev/null --error=/dev/null --time=04:30:00 --mem=2g --array=1-$N_SIMS --job-name=results_gen_tv_$SUBDIR_NAME --wait R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME $LAMBDA_COMPONENT" results_gen_lambda_tv.R results_gen_tv_$SUBDIR_NAME.Rout &&

sbatch --output=/dev/null --error=/dev/null --job-name=analyse_$SUBDIR_NAME R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME $LAMBDA_COMPONENT" analyse_results_lambda.R analyse_results_$SUBDIR_NAME.Rout
