from chunk import chunker
import cyvcf2
import sys
import json

def contig_bounds(vcf):
    def getbounds(vcf, cname):
       contig = list(vcf(cname))
       bounds = [contig[x].POS for x in [0, -1]]
       return bounds
    cnames = vcf.seqnames
    return {x: getbounds(vcf, x) for x in cnames}

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
