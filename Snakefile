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

# Qualimap
annot_gtf = config['gtf']
qualimap_container = config['qualimap_container']
qualimap_mem = config['memory_size']

# MultiQC
multiqc_container = config['multiqc_container']

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

def getQualimapProtocol(wildcards):
    strandinfo = samples.loc[(wildcards.sample_id), 'strandedness']
    if strandinfo == 'reverse':
        return 'strand-specific-reverse'
    elif strandinfo == 'forward':
        return 'strand-specific-forward'
    else:
        return 'non-strand-specific'

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
        expand("rsem/{sample_id}.genes.results", sample_id = samples['id']),
        expand("qc/qualimap/{sample_id}", sample_id = samples['id']),
        "qc/multiqc_report.html",
        "qc/validateFastq_summary.csv"

rule star:
    input:
        unpack(get_fqs),
        genomedir=star_index
    output:
        outdir=directory("star/{sample_id}/"),
        outgenbam="star/{sample_id}/Aligned.sortedByCoord.out.bam",
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
        "rsem/{sample_id}.genes.results",
        "rsem/{sample_id}.isoforms.results",
        "rsem/{sample_id}.transcript.bam"
    threads: rsem_threads
    singularity: rsem_container
    params: 
        ref=rsem_ref,
        outprefix="rsem/{sample_id}",
        strand=lambda wildcards: getRSEMStranding(wildcards)
    shell:
        """
        rsem-calculate-expression --paired-end -p {threads} \
            --alignments {params.strand} {input} {params.ref} {params.outprefix}
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
        biopet-validatefastq {params} | tee {output}
        """

rule checkFastqValidation:
    input:
        expand("qc/validateFastq/{sample_id}.txt", sample_id = samples['id'])
    output:
        "qc/validateFastq_summary.csv"
    script:
        "scripts/summarize_fastqValidation.py"

rule qualimap:
    input:
        bam="star/{sample_id}/Aligned.sortedByCoord.out.bam",
        gtf=annot_gtf
    output:
        directory("qc/qualimap/{sample_id}")
    params:
        mem=qualimap_mem,
        paired=lambda wildcards, input: f"-pe" if isPE(wildcards) else "",
        protocol=lambda wildcards: getQualimapProtocol(wildcards)
    singularity: qualimap_container
    shell:
        """
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} \
            -outdir {output} --java-mem-size={params.mem} \
            -s {params.paired} -p {params.protocol}
        """

rule multiqc:
    input:
        star="star",
        rsem="rsem",
        qc="qc/qualimap"
    output:
        "qc/multiqc_report.html",
        directory("qc/multiqc_data")
    params:
        outdir = "qc/"
    singularity: multiqc_container
    shell:
        """
        multiqc {input} -o {params.outdir}
        """
