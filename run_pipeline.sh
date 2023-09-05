#!/bin/bash
#SBATCH --job-name=star-rsem
#SBATCH --output=/scratch/groups/carilee/nmsc-star/log
#SBATCH --nodes=1
#SBATCH --time=02-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=250
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu

set -e
cd $GROUP_SCRATCH/nmsc-star
echo "Starting snakemake..."
snakemake --cluster-config cluster.json -j 499 \
    --rerun-incomplete \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'

