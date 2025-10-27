#!/bin/bash
#SBATCH --job-name=qupath
#SBATCH --partition=gpu4_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --gres=gpu:2                  
#SBATCH --time=8:00:00               
#SBATCH --output=qupath_%j.out
#SBATCH --error=qupath_%j.err

# Load required modules
module load qupath
module load cuda91/toolkit/9.1.85
nvidia-smi

# Set Java memory options for QuPath
export JAVA_OPTS="-Xmx200g -XX:MaxDirectMemorySize=32g"
export _JAVA_OPTIONS="-Xmx200g"

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "GPUs: $CUDA_VISIBLE_DEVICES"
echo "Start time: $(date)"

# Define file paths
CZI_FILE="/gpfs/data/zwicklab/raz/Raz/07282025_MolecularInstruments_Grp1_compressed.czi"
OUTPUT_DIR="/gpfs/data/zwicklab/raz/Raz/qupath_output_${SLURM_JOB_ID}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run QuPath in GUI mode (requires X11 forwarding)
/gpfs/share/apps/qupath/0.4.3/raw/QuPath/bin/QuPath "$CZI_FILE" &
wait