#!/bin/bash
#SBATCH --partition batch
#SBATCH --qos batch_default
#SBATCH --job-name=correxpv2    
#SBATCH --output=log_%A_%a_%x.out      # Output file, with job name, job ID, and array task ID
#SBATCH --error=log_%A_%a_%x.err       # Error file
#SBATCH --mail-type BEGIN,FAIL,END
#SBATCH --mail-user jcbacong@up.edu.ph
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


# # Define array of datasets
# dataset=("genex_sepsis_only.ap_mi.csv" "genex_shock_only.ap_mi.csv" "genex_diseased.ap_mi.csv" "genex_healthy.ap_mi.csv")
# powers=(10 10 5 5)

TASK_ID=$SLURM_ARRAY_TASK_ID

# # Check if TASK_ID is within bounds of dataset array
# if [ "$TASK_ID" -ge "${#dataset[@]}" ]; then
#    echo "Task ID $TASK_ID exceeds available dataset range."
#    exit 1
# fi

# Get the base name for the task

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

moduleAList=("modules_healthy_scalefree10.csv" "modules_nonshock_scalefree5.csv")
moduleBList=("modules_diseased_scalefree6.csv" "modules_shock_scalefree8.csv")
exprPath="genexdata/"
adjPath="adjmats/"
exprAList=("genex2_nonsepsis.txt" "genex2_preshock.txt")
exprBList=("genex2_sepsis.txt" "genex2_shock.txt")
adjBList=("dissTOM_diseased_dcor6.txt" "dissTOM_shock_dcor8.txt")
adjAList=("dissTOM_healthy_dcor10.txt" "dissTOM_nonshock_dcor5.txt")
outList=("healthyVsdiseased-corrMat-serial" "nonshockVsshock-corrMat-serial")

moduleA=${moduleAList[$TASK_ID]}
moduleB=${moduleBList[$TASK_ID]}
exprA=${exprAList[$TASK_ID]}
exprB=${exprBList[$TASK_ID]}
adjA=${adjAList[$TASK_ID]}
adjB=${adjBList[$TASK_ID]}
output=${outList[$TASK_ID]}

echo "[Task ID: $TASK_ID] Running python script for $moduleA, $moduleB, $exprPath$exprA, $exprPath$exprB ..."

# MCP calculation
start_time=$(date +%s.%N)
python -u correxp-v2.py --moduleA $moduleA --moduleB $moduleB --exprA ${exprPath}${exprA} --exprB ${exprPath}${exprB} --adjA ${adjPath}${adjA} --adjB ${adjPath}${adjB} --output $output
end_time=$(date +%s.%N)

run_time=$(echo "$end_time - $start_time" | bc)
echo "[Task ID: $TASK_ID] Finished histogram calculation $(date)"
echo "[Task ID: $TASK_ID] Total runtime (sec): $run_time"

