#!/bin/bash
#SBATCH --job-name=star-rsem
#SBATCH --output=/scratch/users/tbencomo/nmsc-star/log2
#SBATCH --nodes=1
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu

# Run in non-slurm mode to make validateFastq work

set -e
cd $(pwd)
echo "Starting snakemake..."
snakemake --use-singularity -j 24

