#!/bin/bash
#SBATCH --mem=120g
#SBATCH --cpus-per-task=80
#SBATCH --time=3:00:00
#SBATCH --gres=lscratch:120
#SBATCH --job-name="testing2_cross-study_cluster_matching"

SHEET_PATH=""
TMPDIR=/lscratch/$SLURM_JOB_ID
RESULTS_PATH="../results/"

# read args from CL
while [[ $# -gt 0 ]]; do
    case $1 in 
        --sheet_path) SHEET_PATH="$2"; shift 2;;
        *) echo "Warning: Unknown parameter passed: $1"; shift 1;;
    esac
done

if [[ -z "$SHEET_PATH" ]]; then
    echo "Error: --sheet_path is required."
    exit 1
fi

# activate conda env
source myconda; conda activate frmatch_env

mkdir -p "$RESULTS_PATH"
# INCLUDE RUN LOG WRITING LOGIC
RUN_LOG="${RESULTS_PATH}/run_config_${SLURM_JOB_ID}.txt"
echo "=== Pipeline Run Configuration ===" > "$RUN_LOG"
echo "Date: $(date)" >> "$RUN_LOG"
echo "Slurm Job ID: $SLURM_JOB_ID" >> "$RUN_LOG"
echo "Allocated Cores: $SLURM_CPUS_PER_TASK" >> "$RUN_LOG"
echo "Allocated Memory: ${SLURM_MEM_PER_NODE}MB" >> "$RUN_LOG"
echo "----------------------------------" >> "$RUN_LOG"
echo "Data ID: $DATA_ID" >> "$RUN_LOG"
echo "Sample Sheet Path: $SHEET_PATH" >> "$RUN_LOG"
echo "==================================" >> "$RUN_LOG"

#### DATA INGEST AND PIPELINE ####
python scripts/ingest.py --sheet_path "$SHEET_PATH" --tmpdir "$TMPDIR"

#################### FOR TESTING PURPOSES ####################
cd $TMPDIR || exit

subdirs=(*/)
first_subdir="${subdirs[0]}"

cd "$first_subdir" || exit
pwd

for file in *; do
    if [ -f "$file" ]; then
        FILE_SIZE=$(ls -lh "$file" | awk '{print $5}')
        echo "--- File: $file, Size: $FILE_SIZE ---"
        head -n 5 "$file"
        echo -e "\n"
    fi
done
#################### FOR TESTING PURPOSES ####################

ROOT_DIR=$(dirname "$SLURM_SUBMIT_DIR")
Rscript scripts/run_FR-Match.R "$ROOT_DIR" "$TMPDIR" "$RESULTS_PATH" "$SLURM_CPUS_PER_TASK"