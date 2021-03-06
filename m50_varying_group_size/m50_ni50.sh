#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=8g
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=srosin@live.unc.edu
module load r/4.0.1  

#define variables
RAN_SEED=2031
M=50
N_I=50
N_SIMS=1000
SUBDIR_NAME=M50_NI50

srun --output=/dev/null --error=/dev/null R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME" dataset_gen.R dataset_gen_$SUBDIR_NAME.Rout

sbatch --output=/dev/null --error=/dev/null --time=8:00:00 --mem=2g --array=1-$N_SIMS --job-name=results_gen_$SUBDIR_NAME --wait R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME" results_gen.R results_gen_$SUBDIR_NAME.Rout &&

sbatch --output=/dev/null --error=/dev/null --job-name=analyse_$SUBDIR_NAME R CMD BATCH --no-save --no-restore "--args $RAN_SEED $M $N_I $N_SIMS $SUBDIR_NAME" analyse_results.R analyse_results_$SUBDIR_NAME.Rout
