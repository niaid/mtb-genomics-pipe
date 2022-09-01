#!/bin/bash

# submit this file with:  qsub submit.sh

# run from current working directory
#$ -cwd
#$ -m be

# log dirs
#$ -e ./log/submit_log/
#$ -o ./log/submit_log/

snake_log=$PWD/log/snake_log
mkdir -p $snake_log

# create qsub command
sbcmd="qsub -l {cluster.mem} -pe threaded {cluster.threads} -e $snake_log -o $snake_log "

module load snakemake/6.0.5-Python-3.9.2 || exit 1

# submit to cluster
snakemake --use-conda --jobs 10 --rerun-incomplete --keep-going --cluster-config cluster.yaml --cluster "$sbcmd" --latency-wait 120 all

# # if directory is locked run:
# snakemake --unlock