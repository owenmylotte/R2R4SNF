#!/bin/bash
#SBATCH --job-name=htgr_fcm_le_nosod
#SBATCH --output=htgr_fcm_le_nosod.out
#SBATCH --nice=19  # Set a high niceness value
#SBATCH --nodes=1  # Request one node
#SBATCH --exclusive=user  # Ensure exclusive use of the node
#SBATCH --cpus-per-task=1  # Request one CPU core per task

module load openmc

# Define the directory where the variation files are stored
run_file_dir="/net/clusterhn.cluster.com/home/omylotte/22-251/"

run_file_path="R2R4SNF/inputs/openmc/HTGR_FCM_LE_NOSOD/HTGR_FCM_LE_NOSOD.py"

# Run casmo4e for the current variation file
python -i "${run_file_dir}${run_file_path}"