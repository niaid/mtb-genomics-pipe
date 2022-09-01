
rule get_genome:
    output:
        f"data/reference/{accession}.fasta"
    log:
        "log/reference/get-genome.log",
    params:
        accession=accession
    conda:
        "../envs/reference.yaml"
    shell:"""
        python scripts/fetch_genome.py -i {params.accession} -o {output}
    """

checkpoint genome_faidx:
    input:
        rules.get_genome.output
    output:
        f"data/reference/{accession}.fasta.fai"
    log:
        "log/reference/genome-faidx.log"
    conda:
        "../envs/reference.yaml"
    shell:"""
        samtools faidx {input}
    """

rule genome_dict:
    input:
        rules.get_genome.output
    output:
        f"data/reference/{accession}.dict"
    log:
        "log/reference/create_dict.log"
    conda:
        "../envs/reference.yaml"
    shell:"""
        samtools dict {input} > {output} 2> {log}
    """

rule bwa_index:
    input:
        rules.get_genome.output
    output:
        multiext(f"data/reference/{accession}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "log/reference/bwa_index.log"
    conda:
        "../envs/reference.yaml"
    shell:"""
        bwa index {input}
    """