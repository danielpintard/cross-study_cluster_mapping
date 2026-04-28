#!/bin/bash
#SBATCH --mem=300g
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:300
#SBATCH --job-name="cross-study_cluster_matching"

SHEET_PATH=""
TMPDIR=/lscratch/$SLURM_JOB_ID
# RESULTS_PATH="../results/${DATA_ID}_results"

# read args from CL
while [[ $# -gt 0 ]]; do
    case $1 in 
        --sheet_path) SHEET_PATH="$2"; shift 2;;
    esac
done

if [[ -z "$SHEET_PATH" ]]; then
    echo "Error: --sheet_path is required."
    exit 1
fi

# activate conda env
source myconda; conda activate frmatch_env

# mkdir -p "$RESULTS_PATH"

# INCLUDE RUN LOG WRITING LOGIC

#### DATA INGEST ####
python scripts/ingest.py --sheet_path "$SHEET_PATH" --tmpdir "$TMP_DIR"

