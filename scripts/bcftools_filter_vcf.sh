#!/bin/bash

# submit this file with:  qsub bcftools_filter_vcf.sh

# run from current working directory
#$ -cwd

#$ -m be
#$ -M brendan.jeffrey@nih.gov

# log dirs
#$ -e ./log
#$ -o ./log
# mkdir -p ./log

#$ -l h_vmem=24G

# paths
vcf_dir='../output/pilon/merged'
ppe_coords='../data/external/H37Rv_NC_000962_PEPPE_coords.bed'
results_dir='../output/pilon/merged/georgia'

sample="georgia1296_2022-03-31"
# sample_out='georgia_lin2_clade2_2022-05-04'


module load bcftools

## remove PPE genes by coordinates, only want snps
# bcftools view -v snps -T ^${ppe_coords} ${vcf_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${vcf_dir}/${sample}.noPEPPE.snps.norm.vcf.gz -O z
# bcftools view -v snps -T ^${ppe_coords} ${vcf_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${vcf_dir}/${sample}.noPEPPE.snps.norm.vcf

# ## to REMOVE samples from a VCF
# bcftools view -S ^../georgia_remove_2022-04-25.txt -c 1 ${vcf_dir}/${sample}.noPEPPE.snps.norm.vcf.gz -o ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.vcf

# ## to INCLUDE samples from VCF
# bcftools view -S ../georgia_lin2_clade2_2022-05-04.txt -c 1 ${vcf_dir}/${sample}.noPEPPE.snps.norm.vcf.gz -o ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.vcf
# # bcftools view -S ../georgia_lin2_2022-04-28.txt -c 1 ${vcf_dir}/${sample}.noPEPPE.snps.norm.vcf.gz -o ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.vcf.gz -O z

# ### remove those sites where all samples have a variant - used when removing reference to make multifasta - value = 2 * number of samples
# bcftools view -i 'INFO/AC < 280' ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.vcf -o ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.linspecific.vcf

# ## SNP panel
# module purge
# module load anaconda3
# python multiVCF_to_multifasta_v0.4.py -i ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.linspecific.vcf -o ${results_dir}/${sample_out}.noPEPPE.snps.fa

# # # clean
# rm -f ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.vcf ${vcf_dir}/${sample_out}.noPEPPE.dB.filt.snps.norm.TEMP.linspecific.vcf


# ## vcf with lineage specific samples, include PEPPE, all variants SNPs and indels

# sample_out='georgia_lin2_2022-05-04'
# bcftools view -S ../georgia_lin2_2022-04-28.txt -c 1 ${vcf_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${vcf_dir}/${sample_out}.dB.filt.vars.norm.vcf.gz -O z

# sample_out='georgia_lin2_clade1_2022-05-04'
# bcftools view -S ../georgia_lin2_clade1_2022-05-04.txt -c 1 ${vcf_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${vcf_dir}/${sample_out}.dB.filt.vars.norm.vcf.gz -O z

# sample_out='georgia_lin2_clade2_2022-05-04'
# bcftools view -S ../georgia_lin2_clade2_2022-05-04.txt -c 1 ${vcf_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${vcf_dir}/${sample_out}.dB.filt.vars.norm.vcf.gz -O z

## vcf with lineage specific samples, REMOVE PEPPE, all variants SNPs and indels
bcftools view -v snps,indels -T ^${ppe_coords} ${results_dir}/${sample}.dB.filt.vars.norm.vcf.gz -o ${results_dir}/${sample}.noPEPPE.vars.norm.vcf.gz -O z

sample_out='georgia_lin2_2022-05-04'
bcftools view -S ../georgia_lin2_2022-04-28.txt -c 1 ${results_dir}/${sample}.noPEPPE.vars.norm.vcf.gz -o ${results_dir}/${sample_out}.noPEPPE.dB.filt.vars.norm.vcf.gz -O z

sample_out='georgia_lin2_clade1_2022-05-04'
bcftools view -S ../georgia_lin2_clade1_2022-05-04.txt -c 1 ${results_dir}/${sample}.noPEPPE.vars.norm.vcf.gz -o ${results_dir}/${sample_out}.noPEPPE.dB.filt.vars.norm.vcf.gz -O z

sample_out='georgia_lin2_clade2_2022-05-04'
bcftools view -S ../georgia_lin2_clade2_2022-05-04.txt -c 1 ${results_dir}/${sample}.noPEPPE.vars.norm.vcf.gz -o ${results_dir}/${sample_out}.noPEPPE.dB.filt.vars.norm.vcf.gz -O z



