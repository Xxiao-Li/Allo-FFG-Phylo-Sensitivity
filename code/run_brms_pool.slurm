#!/bin/bash
#SBATCH --job-name=pooled_m1m8
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=3-4
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=xiaoxiao.li@sheffield.ac.uk


module load R/4.4.1-foss-2022b
mkdir -p log_pooled_m1m8


export TASK_ID=${SLURM_ARRAY_TASK_ID}       
export CHAINS=4
export CORES=4
export THREADS=4
export GRAINSIZE=0

export ITER=2000
export WARMUP=1000
export ADAPT_DELTA=0.99
export MAX_TREEDEPTH=18
export SEED=1007

export RESP=log_rel_sensi

Rscript run_pooled_taxa_models.R > rout/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout 2>&1
