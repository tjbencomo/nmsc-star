#!/bin/bash
#SBATCH --job-name=star-rsem
#SBATCH --output=/scratch/users/tbencomo/nmsc-star/log
#SBATCH --nodes=1
#SBATCH --time=01-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=500
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu

set -e
cd $(pwd)
echo "Starting snakemake..."
snakemake --cluster-config cluster.json -j 499 \
    --rerun-incomplete \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'

