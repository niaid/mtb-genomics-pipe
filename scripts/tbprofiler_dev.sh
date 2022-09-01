#!/bin/bash

# submit this file with:  qsub submit.sh

# run from current working directory
#$ -cwd

#$ -m be
#$ -M brendan.jeffrey@nih.gov

# log dirs
#$ -e ./log
#$ -o ./log

reads='/hpcdata/bcbb/jeffreybm/projects/tb_portal/TBportals_genomics_v2_dev/data/raw'
output='/hpcdata/bcbb/jeffreybm/projects/tb_portal/TBportals_genomics_v2_dev/test/tbprofiler'

sample='SRR16087946'

module load tbprofiler/4.1.1_python3
tb-profiler profile --threads 8 \
    --caller pilon \
    --read1 ${reads}/${sample}_1.fastq.gz --read2 ${reads}/${sample}_2.fastq.gz \
    --prefix ${sample} \
    --dir ${output} --txt

# tb-profiler profile --threads 8 \
#     --caller pilon \
#     --read1 ${sample}_1.fastq.gz --read2 ${sample}_2.fastq.gz \
#     --prefix ${sample} \
#     --dir ${output}--txt 2> {log}