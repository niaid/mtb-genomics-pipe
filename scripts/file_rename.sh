#!/bin/bash

# submit this file with:  qsub submit.sh

# run from current working directory
#$ -cwd

#$ -m be
#$ -M brendan.jeffrey@nih.gov

# log dirs
#$ -e ./log
#$ -o ./log
# mkdir -p ./log


# paths
read_dir="../data/raw_temp"

# rename samples
rename .filtered_2nd.1.fastq.gz _1.fastq.gz ${read_dir}/*.filtered_2nd.1.fastq.gz
rename .filtered_2nd.2.fastq.gz _2.fastq.gz ${read_dir}/*.filtered_2nd.2.fastq.gz