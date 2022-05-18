import pandas as pd
from pathlib import Path

configfile: "config.yaml"
samples_fp = config['samples']
samples = pd.read_csv(samples_fp, dtype=str)

# Logs
slurm_logdir = config['slurm_log_dir']
logpath = Path(slurm_logdir)
logpath.mkdir(parents=True, exist_ok=True) 

# STAR
star_index = config['star_index']
star_threads = config['star_threads']
star_container = config['star_container']

# RSEM
rsem_ref = config['rsem_ref']
rsem_threads = config['rsem_threads']
rsem_container = config['rsem_container']

# Set sample indices
samples['id'] = samples['patient'] + '-' + samples['condition']
samples = samples.set_index(["id"], drop=False)
samples = samples.sort_index()

validateFastq_container = config['validateFastq_container']

star_index = config['star_index']

def get_fqs(wildcards):
    if isPE(wildcards):
        return {
            'fq1' : samples.loc[(wildcards.sample_id), 'fq1'],
            'fq2' : samples.loc[(wildcards.sample_id), 'fq2']
        }
    else:
        return {'fq' : samples.loc[(wildcards.sample_id), 'fq1']}

def isPE(wildcards):
    if pd.isna(samples.loc[wildcards.sample_id, 'fq2']):
        return False
    else:
        return True

def isGzZipped(fq):
    return fq.endswith('.gz')

def getRSEMStranding(wildcards):
    strandedness = samples.loc[(wildcards.sample_id), 'strandedness']
    if strandedness == 'reverse':
        return '--forward-prob 0'
    elif strandedness == 'forward':
        return '--forward-prob 1'
    elif strandedness == 'none':
        return '--forward-prob .5'
    else:
        raise ValueError(f"strandedness value {strandedness} not recognized")

rule targets:
    input:
        expand("rsem/{sample_id}", sample_id = samples['id']),
        "qc/validateFastq_summary.csv"

rule star:
    input:
        unpack(get_fqs),
        genomedir=star_index
    output:
        outdir=directory("star/{sample_id}/"),
        outbam="star/{sample_id}/Aligned.toTranscriptome.out.bam"
    threads: star_threads
    singularity: star_container
    params: 
        readcmd = lambda wildcards, input: '--readFilesCommand zcat' if isGzZipped(input.fq1) else '',
        outdir = "star/{sample_id}/"
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {input.genomedir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.outdir} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            {params.readcmd}
        """

rule rsem:
    input:
        "star/{sample_id}/Aligned.toTranscriptome.out.bam"
    output:
        directory("rsem/{sample_id}/")
    threads: rsem_threads
    singularity: rsem_container
    params: 
        ref=rsem_ref,
        outdir="rsem/{sample_id}/",
        strand=lambda wildcards: getRSEMStranding(wildcards)
    shell:
        """
        mkdir {params.outdir}
        rsem-calculate-expression --paired-end -p {threads} \
            --alignments {params.strand} {input} {params.ref} {params.outdir}
        """

rule validateFastq:
    input:
        unpack(get_fqs)
    output:
        "qc/validateFastq/{sample_id}.txt"
    params:
        fqs = lambda wildcards, input: f"-i {input.fq}" if not isPE(wildcards) else f"-i {input.fq1} -j {input.fq2}"
    singularity: validateFastq_container
    shell:
        """
        biopet-validatefastq {params} 2>&1 | tee {output}
        """

rule checkFastqValidation:
    input:
        expand("qc/validateFastq/{sample_id}.txt", sample_id = samples['id'])
    output:
        "qc/validateFastq_summary.csv"
    script:
        "scripts/summarize_fastqValidation.py"
