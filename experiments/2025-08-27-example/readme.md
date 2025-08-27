# Example experiment

## Directory guide

- `config/` - Contains configuration files for jobs (e.g., `Aagos.cfg`)
- `gen-slurm.py` - Python script that generates the slurm job submission files for this experiment.
- `run-gen-slurm.sh` - Bash script that runs `gen-slurm.py` (because `gen-slurm.py` has a bunch of parameters that can be annoying to type into the commandline; easier to write them out in a script)
- `local-run-gen-slurm.sh` - Bash script that runs `gen-slurm.py`, but configured to run on your local machine for testing.

## To run an experiment

1) Copy the compiled `Aagos` executable into this experiment directory's `hpc/config` folder.
2) Run the `run-gen-slurm.sh` script to generate the slurm submission scripts.

    Note: You will need to customize this script (e.g., edit to configure correct file paths, etc)

3) Submit the generated slurm files (should be created inside of whatever was configured for as the `JOB_DIR`)