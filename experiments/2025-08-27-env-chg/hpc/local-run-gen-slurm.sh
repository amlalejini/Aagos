#!/usr/bin/env bash

REPLICATES=100
EXP_SLUG=2025-08-27-env-chg
SEED_OFFSET=10000
JOB_TIME=8:00:00
JOB_MEM=8G
PROJECT_NAME=Aagos
RUNS_PER_SUBDIR=950
USERNAME=lalejina

REPO_DIR=/Users/lalejina/devo_ws/${PROJECT_NAME}
HOME_EXP_DIR=${REPO_DIR}/experiments/${EXP_SLUG}

DATA_DIR=${HOME_EXP_DIR}/hpc/test/data
JOB_DIR=${HOME_EXP_DIR}/hpc/test/jobs
CONFIG_DIR=${HOME_EXP_DIR}/hpc/config

# (1) Activate appropriate Python virtual environment
source ${REPO_DIR}/pyenv/bin/activate

# (2) Generate slurm script
#   - This will generate an events file for each run
python3 gen-slurm.py \
  --runs_per_subdir ${RUNS_PER_SUBDIR} \
  --time_request ${JOB_TIME} \
  --mem ${JOB_MEM} \
  --data_dir ${DATA_DIR} \
  --config_dir ${CONFIG_DIR} \
  --repo_dir ${REPO_DIR} \
  --replicates ${REPLICATES} \
  --job_dir ${JOB_DIR} \
  --seed_offset ${SEED_OFFSET}