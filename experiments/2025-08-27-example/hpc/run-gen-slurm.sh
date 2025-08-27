#!/usr/bin/env bash

REPLICATES=100
EXP_SLUG=2025-08-27-example
SEED_OFFSET=100000
JOB_TIME=8:00:00
JOB_MEM=8G
PROJECT_NAME=Aagos
RUNS_PER_SUBDIR=950
USERNAME=lalejina # <-- CHANGE THIS TO YOUR OWN ACCOUNT USERNAME!

REPO_DIR=/mnt/home/${USERNAME}/devo_ws/${PROJECT_NAME} # <-- CHANGE THIS to where ever you have this repository stored on your account
REPO_SCRIPTS_DIR=${REPO_DIR}/scripts
HOME_EXP_DIR=${REPO_DIR}/experiments/${EXP_SLUG}

DATA_DIR=/mnt/projects/${USERNAME}_project/${PROJECT_NAME}/${EXP_SLUG}
JOB_DIR=${DATA_DIR}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/hpc/config

# (1) Activate appropriate Python virtual environment
source ${REPO_DIR}/hpc-env/clipper-hpc-env.sh
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