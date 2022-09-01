
# modules
import vcf
from collections import defaultdict
import argparse

"""
python script to parse a snpEff annotated vcf and output a multi fasta of SNPs.  
used in snakemake pipeline

uses a list of PPE locus tags - change in future sxcript to take in a VCF where the PPE SNPs removed in input VCF
"""

# parse command line arguments
parser = argparse.ArgumentParser(description='Extract SNPs from snpEff annotated vcf file and output multi fasta')
parser.add_argument('-d', dest='in_peppe_locus', required=True, help='file containg PE/PPE locus tags')
parser.add_argument('-i', dest='in_vcf', required=True, help='in annotated vcf snp file')
parser.add_argument('-o', dest='out_multi_fasta', required=True, help='')
args = parser.parse_args()

# in PE/PPE locus_tag
in_PEPPE = args.in_peppe_locus

# input files
in_vcf_norm = args.in_vcf

# output files
out_multi_fasta = args.out_multi_fasta


# # in PE/PPE locus_tag
# in_PEPPE = "/Users/jeffreybm/Documents/Genomics/TB_portal/data/reference/H37Rv_PEPPE_locus_tags.txt"

# # input files
# in_vcf_norm = "/Users/jeffreybm/Documents/Genomics/TB_portal/results/variant_calls/variants_Amb/romania_tb_merge.variants.DbN.anno.snps.vcf"

# # output files
# out_multi_fasta = "/Users/jeffreybm/Documents/Genomics/TB_portal/results/variant_calls/romania_tb_merge.variants.snps.multiseq.fa"


# PE/PPE locus_tags to list
with open(in_PEPPE) as f:
    PE_PPE_locus_tags = [locus_tag.rstrip() for locus_tag in f]
    
# iupac codes for ambiguous
iupac_dict = {"AG":"R", "GA":"R",
              "CT":"Y", "TC":"Y",
              "GC":"S", "CG":"S",
              "AT":"W", "TA":"W",
              "GT":"K", "TG":"K",
              "AC":"M", "CA":"M"}

# open vcf instance
vcf_reader = vcf.Reader(open(in_vcf_norm, 'r'))

# samples
my_samples = vcf_reader.samples

# initialize dict holding seqs
sample_seqs_dict = defaultdict(list)

# dict to hold number of ambiguous calls per sample
sample_amb_count = {}

# dict to hold number of total variant sites per sample
sample_variant_count = {}

counter = 0
record_count = 1
snp_count = 0

# parse vcf
for record in vcf_reader:   
    if len(record.INFO) > 1:
        
        ### 
        record_count += 1
        ### 
        
        # temporary dictionary that will hold all calls for this site
        temp_call_dict = {}
        if len(record.REF) == 1:
            
            # only need to deal with the first annotation field
            annotation = record.INFO["ANN"][0]

            # parse annotation fields, replace empty info with NA
            anno_fields = ["NA" if v is "" else v for v in annotation.split("|")]

            # extract locus_tag
            s_locus_tag = anno_fields[4]

            # ignore PE/PPE genes
            if s_locus_tag not in PE_PPE_locus_tags:

                # parse each sample, no_call_count to make sure no samples have 'no_call'
                no_call_count = 0
                total_variant_count = 0
                ambiguous_count = 0
                for sample in my_samples:
                    s_rec = record.genotype(sample)

                    # if there is no call for this sample
                    if not s_rec.called:
                        no_call_count += 1
                        
                    else:
                        # increment total variant counter for each sample
                        if s_rec.is_variant:
                            total_variant_count += 1
                            sample_variant_count[sample] = sample_variant_count.get(sample, 0) + 1

                        # if ambiguous because of mixed population
                        if s_rec.is_het:
                            ambiguous_count += 1
                            pot_bases = s_rec.gt_bases.replace("/","")
                            s_call = iupac_dict[pot_bases]

                            # increment ambiguous call count
                            sample_amb_count[sample] = sample_amb_count.get(sample, 0) + 1

                        # extract base
                        else:
                            s_call = s_rec.gt_bases.split("/")[0]

                        # add call to the temp dict
                        temp_call_dict[sample] = s_call

                # append to sequence list dictionary, if no call in any sample
                homzygous_variant_count = total_variant_count - ambiguous_count
                if no_call_count == 0 and homzygous_variant_count > 0:
                    # increment snp counter for map file
                    snp_count += 1
                    
                    # append each to SNP panel
                    for sample in my_samples:
                        sample_seqs_dict[sample].append(temp_call_dict[sample])
                    
                    # write map file with SNP panel postion to ref genome position
                    #out_snp2record_map.write("{}\t{}\t{}\n".format(snp_count, record.POS))

                    # add reference sequence
                    sample_seqs_dict["MTB_H37Rv"].append(record.REF)
                    
# join list of seqs in dictionary into string for each sample, write to multi_fasta
with open(out_multi_fasta, "w") as f:
    for sample in my_samples:
        s_sequence = "".join(sample_seqs_dict[sample])
        f.write(">{}\n{}\n".format(sample, s_sequence))
    
    #add reference genome
    ref_name = "MTB_H37Rv"
    ref_sequence = "".join(sample_seqs_dict[ref_name])
    f.write(">{}\n{}\n".format(ref_name, ref_sequence))
    


