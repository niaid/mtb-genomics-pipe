
rule get_fastq_pe:
    """
    Retrieve paired-read FASTQ files from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw/{sample}_1.fastq",
        "data/raw/{sample}_2.fastq"
    params: sra_dir="data/raw"
    log:
        "log/SRA/{sample}.log"
    threads: 4
    conda:
        "../envs/process_reads.yaml"
    shell:"""
        fasterq-dump -t /hpcdata/scratch --threads {threads} {wildcards.sample} --outdir {params.sra_dir} 2> {log}
    """

rule compress_reads:
    input:  
        reads=get_fastq
    output: 
        out_reads=["data/raw/{sample}_1.fastq.gz", "data/raw/{sample}_2.fastq.gz"]
    shell:"""
        gzip {input} 
    """

rule trim_reads_pe:
    input:
        unpack(get_fastq_gz)
    output:
        r1="data/trim/{sample}.trim.1.fastq.gz",
        r2="data/trim/{sample}.trim.2.fastq.gz"
    params:
        trimmer=" ".join(config["params"]["trimmomatic"]["trimmer"]),
        adapter=config["params"]["trimmomatic"]["illumina_clip"],
        r1_unpaired="data/trim/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="data/trim/{sample}.2.unpaired.fastq.gz"
    log:
        "log/trimmomatic/{sample}.log"
    threads: 8
    conda:
        "../envs/process_reads.yaml"
    shell:"""
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {params.r1_unpaired} \
            {output.r2} {params.r2_unpaired} \
            ILLUMINACLIP:{params.adapter} {params.trimmer} 2>&1 | tee {log}

        rm -f {params.r1_unpaired} {params.r2_unpaired}
    """

    