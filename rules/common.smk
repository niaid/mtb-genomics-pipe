shell.prefix("set -eo pipefail; ")

##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),
    norm_level="norm|norm.decom"

##### Helper functions #####
def get_fastq(wildcards):
    return expand("data/raw/{sample}_{group}.fastq", group=[1, 2], **wildcards)

def get_compressed_fastq(wildcards):
    return expand("data/raw/{sample}_{group}.fastq.gz", group=[1, 2], **wildcards)

def get_fastq_gz(wildcards):
    fastqs = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": f"{fastqs.fq1}.gz", "r2": f"{fastqs.fq2}.gz"}
    return {"r1": f"{fastqs.fq1}.gz"}

def get_trimmed_reads(wildcards):
    return expand("data/trim/{sample}.trim.{group}.fastq.gz", group=[1, 2], **wildcards)

def get_trimmed_reads_dict(wildcards):
    # get reads as a dict
    return {"r1": expand("data/trim/{sample}.trim.1.fastq.gz", sample=wildcards.sample),
            "r2": expand("data/trim/{sample}.trim.2.fastq.gz", sample=wildcards.sample)}

def get_sample_bams(wildcards):
    return expand("output/bwa/{sample}.sort.dedup.bam", sample=wildcards.sample)

def get_sample_bai(wildcards):
    return expand("output/bwa/{sample}.sort.dedup.bam.bai", sample=wildcards.sample)

def get_read_group(wildcards):
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:Illumina'".format(sample=wildcards.sample)

def get_filter_info(wildcards):  
    return config["params"]["filtering"]["vcftools"]

def get_all_vcfs(wildcards):
    return expand("output/pilon/individual/{sample}.filt.vcf.gz", sample=samples.index)

def get_all_vcfs_tbi(wildcards):
    return expand("output/pilon/individual/{sample}.filt.vcf.gz.tbi", sample=samples.index)

def get_profiler_res(wildcards):
    return expand("output/tbprofiler_{tbprofiler_caller}/results/{sample}.results.txt", sample=samples.index, tbprofiler_caller=tbprofiler_caller)
    
def get_lorikeet_res(wildcards):
    return expand("output/lorikeet/temp_merge/{sample}.spoligotype", sample=samples.index)
