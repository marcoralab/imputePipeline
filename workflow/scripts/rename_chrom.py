from pyfaidx import Fasta
import json

fasta = snakemake.input['fasta']
bim_in = snakemake.input['bim']
bim_out = snakemake.output['bim']
json_out = snakemake.output['json']

seqs = Fasta(fasta).keys()
if any(x in seqs for x in ['1', 'chrUn_gl000247', 'GL000247.1']):
    mapping = {'23': 'X', '25': 'X', 'XY': 'X', '24': 'Y',
               '26': 'MT', 'M': 'MT', '0M': 'MT', 'chrM': 'MT',
               'chrMT': 'MT'}
    mapping_chr = {'chr{}'.format(x): str(x) for x in range(1, 23)}
    json_ = {'build': 37}
else:
    mapping = {'23': 'chrX', '25': 'chrX', 'XY': 'chrX', '24': 'chrY',
               'X': 'chrX', 'Y': 'chrY', '26': 'chrM', 'M': 'chrM',
               'MT': 'chrM', '0M': 'chrM', 'chrMT': 'chrM'}
    mapping_chr = {str(x): 'chr{}'.format(x) for x in range(1, 23)}
    json_ = {'build': 38,
             'map': {**mapping_chr, 'Y': 'chrY', 'X': 'chrX',
                     'MT': 'chrM', 'M': 'chrM'}}

mapping.update(mapping_chr)

with open(json_out, "w") as f:
    json.dump(json_, f)

with open(bim_in, 'r') as f_in:
    with open(bim_out, 'w') as f_out:
        for line in f_in:
            rec = line.strip().split()
            update = mapping[rec[0]] if rec[0] in mapping else rec[0]
            updated = '\t'.join([update] + rec[1:])
            print(updated, file=f_out)