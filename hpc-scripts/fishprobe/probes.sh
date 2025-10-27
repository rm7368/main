#!/bin/bash
#SBATCH --job-name=probe_design
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=probe_%j.out
#SBATCH --error=probe_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Razvan.Matei@nyulangone.org

# Load conda
module load anaconda3/gpu/new

eval "$(conda shell.bash hook)"

# Activate your environment
conda activate /gpfs/home/rm7368/.conda/envs/probe_design

# Verify environment
echo "Python path: $(which python)"
echo "Conda env: $CONDA_DEFAULT_ENV"
python -c "import sys; print('Python version:', sys.version)"

# Test biopython import
python -c "from Bio import Entrez; print('Biopython import successful')"

# Run the script
python /gpfs/home/rm7368/fishprobe/designer.py --gene-file /gpfs/home/rm7368/fishprobe/genes.txt --email Razvan.Matei@nyulangone.org

echo "Job completed at $(date)"