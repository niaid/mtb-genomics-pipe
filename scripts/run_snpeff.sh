#!/bin/bash

module load bcftools

java -Xmx5G -jar /hpcdata/bcbb/tbportal/genomics/Software/snpEff/snpEff.jar \
    -c /hpcdata/bcbb/tbportal/genomics/Software/snpEff/snpEff.config \
    -noStats -v -no-downstream -no-upstream Mycobacterium_tuberculosis_h37rv all_tb3386_2022-04-13.dB.filt.vars.norm.TEMP.vcf | \
    sed "s/Chromosome/NC_000962.3/g" | bgzip -c > all_tb3386_2022-04-13.dB.filt.vars.norm.ann.vcf.gz