
rule samtools_stats:
    input:
        "mapped/{sample}.sort.dedup.bam"
    output:
        "qc/samtools-stats/{sample}.txt"
    log:
        "log/samtools-stats/{sample}.log"
    shell:"""
        module load samtools
        samtools stats {input} > {output}
    """

rule multiqc:
    input:
        expand(["qc/samtools-stats/{sample}.txt",
                "qc/fastqc/{sample}_fastqc.zip"], sample=samples.index)
    output:
        "qc/multiqc_report.html"
    log:
        "log/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"
