import pandas as pd

configfile: "config.yaml"
samples_fp = config['samples']
samples = pd.read_csv(samples_fp, dtype=str)
star_index = config['star_index']
star_cores = config['star_cores']
star_container = config['star_container']

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

rule targets:
    input:
        expand("star/{sample_id}", sample_id=samples['id'])

rule star:
    input:
        unpack(get_fqs),
        genomedir=star_index
    output:
        directory("star/{sample_id}")
    params:
        nthreads = star_cores * 3
    singularity: star_container
    shell:
        """
        STAR --runThreadN {params.nthreads} \
            --genomeDir {input.genomedir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {output}
        """
        
