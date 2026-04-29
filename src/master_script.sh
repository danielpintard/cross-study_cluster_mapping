#!/bin/bash
#SBATCH --mem=300g
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:300
#SBATCH --job-name="cross-study_cluster_matching"

SHEET_PATH=""
TMPDIR=/lscratch/$SLURM_JOB_ID
RESULTS_PATH="../results/${DATA_ID}_results"

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
python scripts/ingest.py --sheet_path "$SHEET_PATH" --tmpdir "$TMPDIR"
R scripts/run_FR-Match.R "$SLURM_SUBMIT_DIR" "$TMPDIR" "$RESULTS_PATH" 

#################### FOR TESTING PURPOSES ####################
cd $TMPDIR || exit

subdirs=(*/)
first_subdir="${subdirs[0]}"

cd "$first_subdir" || exit
echo pwd

for file in *; do
    if [ -f "$file" ]; then
        echo "--- File: $file ---"
        head -n 5 "$file"
        echo -e "\n"
    fi
done