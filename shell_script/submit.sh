#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100

python3.8 /home/chandro/junk/ROCKSTAR_pid_upid_test/ROCKSTAR_pid_upid_test.py
