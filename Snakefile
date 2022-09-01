import pandas as pd
import numpy as np
from datetime import date
from snakemake.utils import validate
from snakemake.shell import shell
from datetime import date

# function to update samples df
def add_fq(row):
    return [f"data/raw/{row['sample']}_1.fastq", f"data/raw/{row['sample']}_2.fastq"]

###### Config file and sample sheets #####
# config
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# process sample table
samples = pd.read_csv(config["samples"])
samples['fq1'], samples['fq2'] = zip(*samples.apply(lambda row: add_fq(row), axis=1))
samples.set_index("sample", drop=False, inplace=True)

validate(samples, schema="schemas/samples.schema.yaml")

# helper functions
include: "rules/common.smk"

# variables needed
accession = config["ref"]["accession"]
tbprofiler_caller = config["tool_params"]["tbprofiler"]["caller"]

# generate unique output name for merged files
experiment_id = config["experiment_id"]
isodate = config["date"]

merged_sample_name = f"{experiment_id}_{len(samples.index)}_{isodate}"
merged_sample_name_tbprofiler = f"{merged_sample_name}_{tbprofiler_caller}"
norm_level=["norm", "norm.decom"]


##### Target rules #####
multi_vcf = f"results/{merged_sample_name}/{merged_sample_name}.filt.vars.vcf.gz"
tbprofiler_dst = f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.json"
tbprofiler_vars = f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.variants.txt"
tbprofiler_summary = f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.summary.txt"

# temp dev
reference = f"data/reference/{accession}.fasta"
reads = expand("data/trim/{sample}.trim.{read}.fastq.gz", sample=samples.index, read=['1','2'])
ind_vcf = expand("output/pilon/individual/{sample}.filt.vcf.gz", sample=samples.index)

rule all:
    # input: multi_vcf
    input: tbprofiler_dst, tbprofiler_vars, tbprofiler_summary, multi_vcf

##### Modules #####
include: "rules/reference.smk"
include: "rules/process_reads.smk"
include: "rules/mapping.smk"
include: "rules/variants.smk"
include: "rules/tbprofiler.smk"
