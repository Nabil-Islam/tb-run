#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Tue Feb 18 2025 17:18:14 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=sapphire

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="dowload and tbprof"

# The project ID which this job should run under:
#SBATCH --account="punim1703"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

# The amount of memory in megabytes per node:
#SBATCH --mem=53600

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=7-00:00:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

dir="/data/scratch/projects/punim1703/false_negative/tbprof-rerun/workflow/scripts"

# The job command(s):
python3 "$dir/run_tb.py"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

