import json
import argparse

"""
adds top-level data to TBprofiler.json output
Tb-profiler version run date
"""

# cl args
parser = argparse.ArgumentParser(description='Update TBprofiler.json output with TBprofiler version and run date')
parser.add_argument('-i', dest='in_json', required=True, help='merged TBprofiler.json in')
parser.add_argument('-s', dest='in_lorikeet', required=True, help='merged lorikeet in')
parser.add_argument('-o', dest='out_json', required=True, help='updated TBprofiler.json')
parser.add_argument('-d', dest='date', required=True, help='')
parser.add_argument('-v', dest='tbprofiler_version', required=True, help='')
parser.add_argument('-c', dest='tbprofiler_caller', required=True, help='')
args = parser.parse_args()

# gather spoligotypes
spoligo_dict = {}
with open(args.in_lorikeet) as f:
    for line in f:
        line = line.rstrip()
        if not line.startswith("#"):
            fields = line.split("\t")
            sample = fields[0].split(".")[0]
            spoligotype = fields[1]
            spoligo_dict[sample] = spoligotype

# add spoligotypes
data = json.load(open(args.in_json))
for key in data.keys():
        data[key]['spoligotype'] = spoligo_dict[key]

new_top = {}
new_top["tbprofiler_version"] = args.tbprofiler_version
new_top["run_date"] = args.date 
new_top["variant_caller"] = args.tbprofiler_caller 

new_top["tbprofiler_results"] = data

with open(args.out_json, 'w') as f:
    json.dump(new_top, f, indent=2)