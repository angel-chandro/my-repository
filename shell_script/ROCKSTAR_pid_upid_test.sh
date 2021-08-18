#!/bin/sh

#SBATCH --nodes=1                               # How many nodes? In this case only 1 (i.e. brutus)
#SBATCH --account=32cores                       # The queue name we are interested in.
#SBATCH --cpus-per-task=32                      # How many cores per slurm task
#SBATCH --tasks-per-node=1                      # Run a single slurm task, i.e. only one call of sbatch
#SBATCH --time=1-01                             # Time limit in days-hrs:min:sec
#SBATCH --mem=50000                             # Memory per compute node (in MB)
#SBATCH --array=0-4%1                           # Create a job array with 5 tasks (i.e. 0 to 4) 
                                                # where only 1 task is active at the same time
#SBATCH --job-name=mySlurmJob                   # Job name
#SBATCH --output /home/chandro/output/ROCKSTAR_pid_upid_test  # Output filename
#SBATCH --error  /home/chandro/output/ROCKSTAR_pid_upid_test   # Error filename
# Print general info
pwd; hostname; date
#  Run the script
python3.8 junk/ROCKSTAR_pid_upid_test/ROCKSTAR_pid_upid_test.py
date
