#!/usr/bin/env bash
set -eu

JOB_NAME="snakemake_master_process."$(date "+%s")
LOG_DIR="logs"

if [[ ! -d "$LOG_DIR" ]]; then
  echo "Error: Log directory $LOG_DIR does not exist"
  exit 1
fi

MEMORY="12G"
TIME="${TIME:-1w}"
THREADS=2
PROFILE="slurm.punim1703"
DEFAULT_TMP="tmpdir='/tmp'"
JOBS="${JOBS:-150}"   # max SLURM jobs running/queued at once (one sample per job)
CMD="snakemake --profile $PROFILE --jobs $JOBS --default-resources \"$DEFAULT_TMP\" --local-cores $THREADS --use-conda --conda-frontend mamba --keep-going --rerun-incomplete --latency-wait 60 --rerun-triggers mtime"

ssubmit -t "$TIME" -m "$MEMORY" -o "$LOG_DIR"/"$JOB_NAME".o \
  -e "$LOG_DIR"/"$JOB_NAME".e "$JOB_NAME" "$CMD" -- -c "$THREADS" -p sapphire
