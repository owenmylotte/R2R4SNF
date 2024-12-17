#!/bin/bash
#SBATCH --job-name=casmo_variations
#SBATCH --array=1000-1499%12  # Adjust the range according to the number of pin variations
#SBATCH --output=pset6_%a.out
#SBATCH --nice=19  # Set a high niceness value
#SBATCH --nodes=1  # Request one node
#SBATCH --exclusive=user  # Ensure exclusive use of the node
#SBATCH --cpus-per-task=1  # Request one CPU core per task

module load openmc

# Define the directory where the variation files are stored
run_file_dir="/net/clusterhn.cluster.com/home/omylotte/22-251/"

run_file_path="R2R4SNF/inputs/openmc/HTGR_FCM_LE_SOD/HTGR_FCM_LE_SOD.py"

# Run casmo4e for the current variation file
python -i "${run_file_dir}${run_file_path}"