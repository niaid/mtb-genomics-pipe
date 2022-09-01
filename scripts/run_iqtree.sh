#!/bin/bash

# submit this file with:  qsub run_iqtree.sh

# run from current working directory
#$ -cwd

#$ -m be
#$ -M brendan.jeffrey@nih.gov

# log dirs
#$ -e ./log
#$ -o ./log
# mkdir -p ./log

#$ -l h_vmem=48G

# note: need to get the fconst sites
# (genome size[4.4M] - number of SNPs)/6 for the A and T values, G and C are 3X

sample='georgia_lin2_2022-04-28'

aln_dir='/hpcdata/bcbb/jeffreybm/projects/tb_portal/TBportals_genomics_v2_dev/output/multifasta'
iqtree_dir='/hpcdata/bcbb/jeffreybm/projects/tb_portal/TBportals_genomics_v2_dev/output/iqtree'

module load iqtree
# iqtree -s ${aln_dir}/georgia1297_2022-03-28.noPPE.snps.fa -st DNA -m TEST  -bb 1000 -nt AUTO -redo -pre ${iqtree_dir}/georgia1297_2022-03-28.noPPE
# iqtree -s ${aln_dir}/${sample}.noPEPPE.snps.fa -st DNA -m TEST -fconst 766011,1422592,1422592,766011 -bb 1000 -nt AUTO -redo -pre ${iqtree_dir}/${sample}.noPPE

iqtree -s ${aln_dir}/${sample}.noPEPPE.noref.snps.fa -st DNA -m TEST -fconst 168596,1427392,1427392,168596 -bb 1000 -nt AUTO -redo -pre ${iqtree_dir}/${sample}.noPPE.norefs