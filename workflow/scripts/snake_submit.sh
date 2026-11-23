#!/usr/bin/env bash
#SBATCH --job-name=snakemake_master
#SBATCH --output=logs/snmakemake_master_%j.out
#SBATCH --error=logs/snmakemake_master_%j.err
#SBATCH --time=7-00:00:00          # 1 day (format: D-HH:MM:SS)
#SBATCH --cpus-per-task=4          # threads for master process
#SBATCH --mem=12G
#SBATCH --partition=sapphire       # adjust to your HPC partition

# Load modules

# Paths and Singularity binds
BIND_PATH="/data/scratch/projects/punim1703"
SINGULARITY_ARGS="-B $BIND_PATH"

# Resource for downloads
DOWNLOADS=5 # max concurrent downloads at a time
THREADS=4   # threads for the Snakemake master process
DEFAULT_TMP="tmpdir='/tmp'"

# Run Snakemake
snakemake \
  --profile slurm.punim1703 \
  --rerun-incomplete \
  --default-resources "$DEFAULT_TMP" \
  --jobs 200 \
  --local-cores $THREADS \
  --resources downloads=$DOWNLOADS \
  --conda-frontend conda $* \
  --singularity-args "$SINGULARITY_ARGS" \
  "$@"
