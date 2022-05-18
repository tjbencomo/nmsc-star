# star-rsem
Simple Snakemake pipeline to align RNA-Seq reads and quantify counts using STAR and RSEM. 

## Installation
1. Install Snakemake
```
# Using conda
conda install -c biocoonda snakemake
```
2. Download STAR and RSEM indexes or create your own
3. Ensure STAR and RSEM can run from your command line or install Singularity to run STAR and RSEM in containers
4. Specify the pipeline settings in `config.yaml`
5. Specify the samples to process in `samples.csv`

## Sample File
`samples.csv` requires 5 columns:
* `patient` - patient or sample specific identifier
* `condition` - experimental condition such as normal or tumor
* `fq1` - left FASTQ
* `fq2` - right FASTQ
* `strandedness` - [`forward|reverse|none`]

See [this tutorial](https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/) to determine what stranding your FASTQ files use. `forward` matches case A in the tutorial. `reverse` is case B, and `none` is case C. You can also use [`check-strand`](https://github.com/tjbencomo/check-strand) to quickly infer the strand type using Kallisto.

## Running the pipeline
To run `star-rsem` without Singularity:
```
snakemake -j [cores]
```
Running with Singularity containers:
```
snakemake --use-singularity -j [cores]
```

### Cluster Execution
`run_pipeline.sh` is a bash script that executes the pipeline as a SLURM job. 
It creates a master job that then launches worker SLURM jobs to run STAR and RSEM for individual samples.
To use SLURM execution do the following:
1. Modify `run_pipeline.sh` to run on your cluster
2. Edit the `out` and `account` fields  for the default job in `cluster.json`. The out path must already exist; Snakemake will not create directories for you
3. Launch the master job:
```
sbatch run_pipeline.sh $(pwd)
```

## Output
Two directories will be created
* `star/` - contains the output of each STAR run
* `rsem/` - contains the output of each RSEM run
* `qc/` - contains FASTQ quality control checks via `validateFastq`. These are summarized in `validateFastq_summary.csv`

## Troubleshooting
Common issues:
1. STAR needs a lot of RAM, especially for the human genome. Specify the resources in `cluster.json` accordingly
2. Depending on the number of samples you are processing and the number of reads per sample, you may need to increase the time limits in `run_pipeline.sh` and `cluster.json`
