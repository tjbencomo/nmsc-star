#!/bin/bash
#SBATCH --job-name=star-rsem
#SBATCH --output=
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=END
#SBATCH --mail-user=

set -e
cd $1
echo "Starting snakemake..."
snakemake --cluster-config cluster.json -j 499 \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'

