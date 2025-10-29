#!/usr/bin/env bash

REPLICATES=10
EXP_SLUG=2025-10-29-arch-enumeration
SEED_OFFSET=20000
JOB_TIME=00:10:00
JOB_MEM=8G
PROJECT_NAME=Aagos
RUNS_PER_SUBDIR=950
USERNAME=lalejini
# HPC_ENV_FILE=clipper-hpc-env.sh
HPC_ENV_FILE=msu-hpc-env.sh

REPO_DIR=/mnt/home/${USERNAME}/devo_ws/${PROJECT_NAME} # <-- CHANGE THIS to where ever you have this repository stored on your account
REPO_SCRIPTS_DIR=${REPO_DIR}/scripts
HOME_EXP_DIR=${REPO_DIR}/experiments/${EXP_SLUG}

DATA_DIR=/mnt/scratch/${USERNAME}/data/${PROJECT_NAME}/${EXP_SLUG}
JOB_DIR=${DATA_DIR}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/hpc/config
HPC_ENV_FILEPATH=${REPO_DIR}/hpc-env/${HPC_ENV_FILE}

# (1) Activate appropriate Python virtual environment
source ${HPC_ENV_FILEPATH}
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
  --seed_offset ${SEED_OFFSET} \
  --hpc_account ${ACCOUNT} \
  --hpc_env_file ${HPC_ENV_FILEPATH}