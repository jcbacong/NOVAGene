#!/bin/bash
#SBATCH --job-name=dcor
#SBATCH --output=log_%A_%a_%x.out      # Output file, with job name, job ID, and array task ID
#SBATCH --error=log_%A_%a_%x.err       # Error file
#SBATCH --mail-type BEGIN,FAIL,END
#SBATCH --mem 90G

# Load any necessary modules (e.g., for Conda if it’s not automatically available)
# module load anaconda

# Initialize Conda (if not automatically available in the environment)
module purge
module load anaconda/3-2023.07-2

# Step 2: Activate the new environment
conda activate meowth

# Set stack size to unlimited
ulimit -s unlimited

echo "============ COMPUTE SETTINGS ==========="
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_CPUS_ON_NODE="$SLURM_CPUS_ON_NODE
echo "working directory = "$SLURM_SUBMIT_DIR
echo "========================================="
echo ""


# Define array of datasets
dataset=("genex2_nonsepsis.txt" "genex2_sepsis.txt" "genex2_preshock.txt" "genex2_shock.txt")
output=("dcor_nonsepsis.csv" "dcor_sepsis.csv" "dcor_preshock.csv" "dcor_shock.csv")

TASK_ID=$SLURM_ARRAY_TASK_ID

# Check if TASK_ID is within bounds of dataset array
if [ "$TASK_ID" -ge "${#dataset[@]}" ]; then
   echo "Task ID $TASK_ID exceeds available dataset range."
   exit 1
fi

# Get the base name for the task
DIR="genexdata/"
base_name=${dataset[$TASK_ID]}
output_name=${output[$TASK_ID]}

echo ""
echo "============ TASK CONFIG ==========="
echo "TASK ID: $TASK_ID"
echo "TOTAL NODES: $SLURM_NNODES"
echo "CPUs USED: $SLURM_CPUS_PER_TASK" 
echo "CPUs TOTAL: $SLURM_CPUS_ON_NODE"
echo "INPUT_DIR: $MCP_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "========================================="
echo ""

file=${DIR}${base_name}
echo "[Task ID: $TASK_ID] Running python script for $file ..."

# MCP calculation
start_time=$(date +%s.%N)
python -u calc_dcor.py --input $file --output $output_name
end_time=$(date +%s.%N)

run_time=$(echo "$end_time - $start_time" | bc)
echo "[Task ID: $TASK_ID] Finished dcor calculation $(date)"
echo "[Task ID: $TASK_ID] Total runtime (sec): $run_time"

