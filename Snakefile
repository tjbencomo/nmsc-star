import pandas as pd

configfile: "config.yaml"
samples_fp = config['samples']
samples = pd.read_csv(samples_fp, dtype=str)

star_index = config['star_index']
star_threads = config['star_threads']
star_container = config['star_container']

rsem_ref = config['rsem_ref']
rsem_threads = config['rsem_threads']
rsem_container = config['rsem_container']

samples['id'] = samples['patient'] + '-' + samples['condition']
samples = samples.set_index(["id"], drop=False)
samples = samples.sort_index()

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

def isZipped(fq):
    return fq.endswith('.gz')

rule targets:
    input:
        expand("rsem/{sample_id}", sample_id = samples['id'])
        #expand("star/{sample_id}", sample_id=samples['id'])

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
        readcmd = lambda wildcards, input: '--readFilesCommand zcat' if isZipped(input.fq1) else '',
        outdir = "star/{sample_id}/"
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {input.genomedir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.outdir} \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
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
        outdir="rsem/{sample_id}/"
    shell:
        """
        mkdir {params.outdir}
        rsem-calculate-expression --paired-end -p {threads} \
            --alignments {input} {params.ref} {params.outdir}
        """
     
