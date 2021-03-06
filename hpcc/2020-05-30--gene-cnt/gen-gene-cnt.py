'''
Generate slurm job submission script for 2020-05-18 -- environmental change rate experiment.

See 2020-05-18--env-chg-sweep/README.md for more details.
'''

import argparse, os, sys, errno, subprocess, csv

seed_offset = 500000
default_num_replicates = 50
job_time_request = "4:00:00"
job_memory_request = "8G"
job_name = "gene-cnt"

nk_config = {
    "environment_change": [
        "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0",
        "-CHANGE_MAGNITUDE 64 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 128 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 256 -CHANGE_FREQUENCY 1"
    ]
}

gradient_config = {
    "environment_change": [
        "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0",   # 0
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 4",   # 7
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 2",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 1"
    ]
}

shared_config = {
    "BIT_FLIP_PROB": [
        "-BIT_FLIP_PROB 0.003 -PHASE_2_BIT_FLIP_PROB 0.003",
        "-BIT_FLIP_PROB 0.03 -PHASE_2_BIT_FLIP_PROB 0.03"
    ],

    "GENOME_SIZE": [
        "-NUM_BITS 128 -NUM_GENES 8 -MAX_SIZE 1024",
        "-NUM_BITS 128 -NUM_GENES 32 -MAX_SIZE 1024",
        "-NUM_BITS 128 -NUM_GENES 64 -MAX_SIZE 1024"
    ]
}

base_resub_script = \
"""#!/bin/bash
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=<<TIME_REQUEST>>          # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --array=<<ARRAY_ID_RANGE>>
#SBATCH --mem=<<MEMORY_REQUEST>>        # memory required per node - amount of memory (in bytes)
#SBATCH --job-name <<JOB_NAME>>         # you can give your job a name for easier identification (same as -J)
#SBATCH --account=devolab

########## Command Lines to Run ##########

EXEC=Aagos
CONFIG_DIR=<<CONFIG_DIR>>

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.7.0

<<RESUBMISSION_LOGIC>>

mkdir -p ${RUN_DIR}
cd ${RUN_DIR}
cp ${CONFIG_DIR}/Aagos.cfg .
cp ${CONFIG_DIR}/${EXEC} .

./${EXEC} ${RUN_PARAMS} > run.log

rm Aagos.cfg
rm ${EXEC}

"""

base_run_logic = \
"""
if [[ ${SLURM_ARRAY_TASK_ID} -eq <<RESUB_ID>> ]] ; then
    RUN_DIR=<<RUN_DIR>>
    RUN_PARAMS=<<RUN_PARAMS>>
fi
"""



'''
This is functionally equivalent to the mkdir -p [fname] bash command
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def extract_settings(run_config_path):
    content = None
    with open(run_config_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    content = content[1:]
    configs = [l for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return {param[header_lu["parameter"]]:param[header_lu["value"]] for param in configs}

def is_run_complete(path):
    # (1) Does the run directory exist?
    print(f"    Run dir? {os.path.exists(path)}")
    if not os.path.exists(path): return False
    # (2) If the run directory exists, did the run complete?
    #     Is there a run config file?
    run_config_path = os.path.join(path, "output", "run_config.csv")
    print(f"    Run config? {os.path.exists(run_config_path)}")
    if not os.path.exists(run_config_path): return False
    #    The run config file exists, extract parameters.

    run_params = extract_settings(run_config_path)
    final_gen = run_params["TOTAL_GENS"] # We'll look for this generation in the fitness.csv file
    fitness_file_path = os.path.join(path, "output", "fitness.csv")
    print(f"    Fitness file? {os.path.exists(fitness_file_path)}")
    if not os.path.exists(fitness_file_path): return False

    fitness_contents = None
    with open(fitness_file_path, "r") as fp:
        fitness_contents = fp.read().strip().split("\n")
    if len(fitness_contents) == 0: return False

    header = fitness_contents[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    last_line = fitness_contents[-1].split(",")
    print(f"    len(header) == len(last_line)? {len(header) == len(last_line)}")
    if len(header) != len(last_line): return False

    final_fitness_update = last_line[header_lu["update"]]
    print(f"    {final_fitness_update} =?= {final_gen}")
    if final_fitness_update != final_gen: return False
    return True

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--config_dir", type=str, help="Where is the configuration directory for experiment?")
    parser.add_argument("--replicates", type=int, default=default_num_replicates, help="How many replicates should we run of each condition?")
    parser.add_argument("--query_condition_cnt", action="store_true", help="How many conditions?")

    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    num_replicates = args.replicates

    # Compute all combinations of NK fitness model settings and gradient fitness settings
    nk_combos = [f"{chg} {mut} {genome} -GRADIENT_MODEL 0" for chg in nk_config["environment_change"] for mut in shared_config["BIT_FLIP_PROB"] for genome in shared_config["GENOME_SIZE"] ]
    gradient_combos = [f"{chg} {mut} {genome} -GRADIENT_MODEL 1" for chg in gradient_config["environment_change"] for mut in shared_config["BIT_FLIP_PROB"] for genome in shared_config["GENOME_SIZE"] ]

    # Combine
    combos = gradient_combos + nk_combos
    if (args.query_condition_cnt):
        print("Conditions", combos)
        print(f"Number of conditions: {len(combos)}")
        exit(0)

    # Find complete/incomplete runs.
    num_finished = 0
    resubmissions = []
    for condition_id in range(0, len(combos)):
        condition_params = combos[condition_id]
        print(f"Processing condition: {condition_params}")
        # Run N replicates of this condition.
        for i in range(1, num_replicates+1):
            # Compute seed for this replicate.
            seed = seed_offset + (condition_id * num_replicates) + i
            # Generate run parameters, use to name run.
            run_params = condition_params + f" -SEED {seed}"
            run_name = ("RUN " + run_params).replace("-", "_").replace(" ", "_")
            run_dir = os.path.join(data_dir, run_name)
            # (1) Does the run directory exist?
            print(f"  {run_params}")
            run_complete = is_run_complete(run_dir)
            print(f"    finished? {run_complete}")
            num_finished += int(run_complete)
            if not run_complete: resubmissions.append({"run_dir": run_dir, "run_params": run_params})

    print(f"Runs finished: {num_finished}")
    print(f"Resubmissions: {len(resubmissions)}")

    print("Generating resubmission script...")
    if len(resubmissions) == 0: return

    resub_logic = ""
    array_id = 1
    for resub in resubmissions:
        run_params = resub["run_params"]
        run_logic = base_run_logic
        run_logic = run_logic.replace("<<RESUB_ID>>", str(array_id))
        run_logic = run_logic.replace("<<RUN_DIR>>", resub["run_dir"])
        run_logic = run_logic.replace("<<RUN_PARAMS>>", f"'{run_params}'")

        resub_logic += run_logic
        array_id += 1

    script = base_resub_script
    script = script.replace("<<TIME_REQUEST>>", job_time_request)
    script = script.replace("<<ARRAY_ID_RANGE>>", f"1-{len(resubmissions)}")
    script = script.replace("<<MEMORY_REQUEST>>", job_memory_request)
    script = script.replace("<<JOB_NAME>>", job_name)
    script = script.replace("<<CONFIG_DIR>>", config_dir)
    script = script.replace("<<RESUBMISSION_LOGIC>>", resub_logic)

    with open("gene_size_resub.sb", "w") as fp:
        fp.write(script)

if __name__ == "__main__":
    main()