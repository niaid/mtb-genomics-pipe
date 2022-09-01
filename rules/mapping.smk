
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp("output/bwa/{sample}.sort.bam")
    log:
        "log/bwa_mem/{sample}.log"
    params:
        index=f"data/reference/{accession}.fasta",
        extra=get_read_group
    threads: 8
    conda:
        "../envs/mapping.yaml"
    shell:"""
        bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort - > {output} 2> {log}
    """

rule dedup_bam:
    input:
        "output/bwa/{sample}.sort.bam"   
    output:
        "output/bwa/{sample}.sort.dedup.bam"
    conda:
        "../envs/mapping.yaml"
    shell:"""
        samtools rmdup {input} {output}
        """

rule index_bam:
    input:
        "output/bwa/{sample}.sort.dedup.bam"
    output:
        "output/bwa/{sample}.sort.dedup.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:"""
        samtools index {input}
    """
     