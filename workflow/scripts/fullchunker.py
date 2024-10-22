import json
import sys
import re
import cyvcf2
from chunks import chunker

def contig_bounds(vcf):
    def getbounds(vcf, cname):
       contig = [x for x in vcf(cname)]
       if not contig:
           return None
       bounds = [contig[x].POS for x in [0, -1]]
       return bounds
    cnames = [x for x in vcf.seqnames
              if re.match(r'^(chr){0,1}([12]{0,1}\d|X|Y|M|MT)$', x)]
    allbounds = {x: getbounds(vcf, x) for x in cnames}
    return {k: v for k, v in allbounds.items() if v}

def chunk_contigs(bounds):
    chunks = []
    for contig, range in bounds.items():
        chunks += chunker(contig, 20000000, range[1], range[0], raw=True)
    return chunks

def make_chunk_json(fname, chunks):
    keys = ['chrom', 'from', 'through']
    chunkdict = {x: y for x, y in zip(keys, zip(*chunks))}
    with open(fname, 'w') as file:
        json.dump(chunkdict, file)

if __name__ == "__main__":
    in_name = snakemake.input[0]
    out_name = snakemake.output[0]
    vcf = cyvcf2.VCFReader(in_name)
    bounds = contig_bounds(vcf)
    chunks = chunk_contigs(bounds)
    make_chunk_json(out_name, chunks)
