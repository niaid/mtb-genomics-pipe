
# modules
import vcf
from collections import defaultdict
import argparse
import math
import numpy as np

"""
python script to parse a snpEff annotated vcf and output a multi fasta of SNPs.  
used in snakemake pipeline

"""

# version 0.2 - added functionality to deal with reference sequence having N
# version 0.3 - added indel parsin functionality, reduces indels to shortest presence absence representation
# version 0.4 - fixed reference sequence bug

# parse command line arguments
parser = argparse.ArgumentParser(description='Extract SNPs from snpEff annotated vcf file and output multi fasta')
parser.add_argument('-i', dest='in_vcf', required=True, help='in annotated vcf snp file')
parser.add_argument('-o', dest='out_multi_fasta', required=True, help='')
args = parser.parse_args()

# input files
in_vcf_norm = args.in_vcf

# output files
out_multi_fasta = args.out_multi_fasta

# reference name
ref_name = "MH37Rv_NC_000962"


# # in PE/PPE locus_tag
# in_PEPPE = "/Users/jeffreybm/Documents/Genomics/TB_portal/data/reference/H37Rv_PEPPE_locus_tags.txt"

# # input files
# in_vcf_norm = "/Users/jeffreybm/Documents/Genomics/TB_portal/results/variant_calls/variants_Amb/romania_tb_merge.variants.DbN.anno.snps.vcf"

# # output files
# out_multi_fasta = "/Users/jeffreybm/Documents/Genomics/TB_portal/results/variant_calls/romania_tb_merge.variants.snps.multiseq.fa"

    
## Functions ##
def pad_indels(calls):
    
    # trim calls to shortest unique
    calls = [call[0:len(calls) + 1] for call in calls]
    
    # determine if indels
    lengths = [len(call) for call in calls]
    if len(np.unique(lengths)) > 1:  # indel
        longest = max(calls, key = len)
        calls_pad = [item.ljust(len(longest), '-')[1:len(calls)] for item in calls]
        return calls_pad
    
    else:
        return calls
        
# iupac codes for ambiguous
iupac_dict = {"AG":"R", "GA":"R",
              "CT":"Y", "TC":"Y",
              "GC":"S", "CG":"S",
              "AT":"W", "TA":"W",
              "GT":"K", "TG":"K",
              "AC":"M", "CA":"M"}

# open vcf instance
vcf_reader = vcf.Reader(open(in_vcf_norm, 'r'))

# sample names and number
my_samples = vcf_reader.samples
num_samples = len(my_samples)

# cutoff number for no-calls. either calculate percent or use set
CUT_PERCENT = 0.002
no_call_cutoff = math.floor(CUT_PERCENT * num_samples)

# no_call_cutoff = 2

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

    ### 
    record_count += 1
    ###

    # temporary dictionary that will hold all calls for this site
    temp_call_dict = {}
    # if len(record.REF) == 1:  # dont know why I used the length of the REF before

    # parse each sample, no_call_count to make sure no samples have 'no_call'
    no_call_count = 0
    total_variant_count = 0
    ambiguous_count = 0
    
    # generate lists with calls
    all_calls = []
    all_calls.append(record.REF)
    all_calls = all_calls + record.ALT
    all_calls = [str(call) for call in all_calls]
    
    #pad with - for indels
    all_calls_pad = pad_indels(all_calls)
 
    for sample in my_samples:
        s_rec = record.genotype(sample)
        
        # if there is no call for this sample
        if not s_rec.called:
            no_call_count += 1
            s_call = 'N'.ljust(len(all_calls_pad[0]), '-')
            
            # add call to the temp dict
            temp_call_dict[sample] = s_call

        # if ambiguous because of mixed population
        else:
            #print(record, sample)
            if s_rec.is_variant:
                total_variant_count += 1
                sample_variant_count[sample] = sample_variant_count.get(sample, 0) + 1
                
            if s_rec.is_het:
                
                ambiguous_count += 1
                    
                # call reference if this is an indel site
                if record.is_indel:
                    s_call = all_calls_pad[0]
                    
                # rare case of N in reference and an ambiguous call in sample, call ref  
                elif "N" in all_calls_pad[0]:
                    s_call = all_calls_pad[0]

                else:
                    pot_bases = s_rec.gt_bases.replace("/","")
                    s_call = iupac_dict[pot_bases]

                # increment ambiguous call count
                sample_amb_count[sample] = sample_amb_count.get(sample, 0) + 1

            # extract base
            else:
                genotype = int(s_rec['GT'].split('/')[0])
                s_call = all_calls_pad[genotype]

            # add call to the temp dict
            temp_call_dict[sample] = s_call

    # append to sequence list dictionary, if no call in any sample
    homozygous_variant_count = total_variant_count - ambiguous_count

    if no_call_count <= no_call_cutoff and homozygous_variant_count > 0:
        # increment snp counter for map file
        snp_count += 1

        # append each to SNP panel
        for sample in my_samples:
            sample_seqs_dict[sample].append(temp_call_dict[sample])
            
        # add reference sequence
        sample_seqs_dict[ref_name].append(all_calls_pad[0])     
        
# join list of seqs in dictionary into string for each sample, write to multi_fasta
with open(out_multi_fasta, "w") as f:
    for sample in my_samples:
        s_sequence = "".join(sample_seqs_dict[sample])
        f.write(">{}\n{}\n".format(sample, s_sequence))

    #add reference genome
    ref_sequence = "".join(sample_seqs_dict[ref_name])
    f.write(">{}\n{}\n".format(ref_name, ref_sequence))
    


