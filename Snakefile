'''Snakefile for MIS preparation
   Version 0.8'''

from scripts.parse_config import parser
import socket
import sys

configfile: "config/config.yaml"
shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]

CHROM, SAMPLE, INPATH, keep_command = parser(config)

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

if isMinerva:
    anacondapath = sys.exec_prefix + "/bin"
    shell.prefix(". ~/.bashrc; PATH={}:$PATH; ".format(anacondapath))

if isMinerva:
    com = {'plink': 'plink --keep-allele-order',
           'bcftools': 'bcftools', 'R': 'Rscript',
           'vcftools': 'vcftools'}
    loads = {'plink': 'module load {plink}'.format(plink=config['modules']['plink']),
             'bcftools': 'module load {bcftools}'.format(bcftools=config['modules']['bcftools']),
             'R': 'module load {r}'.format(r=config['modules']['R']),
             'vcftools': 'module load {vcftools}'.format(vcftools=config['modules']['vcftools'])}
else:
    com = {'plink': 'plink --keep-allele-order',
           'bcftools': 'bcftools', 'R': 'Rscript',
           'vcftools': 'vcftools'}
    loads = {'plink': 'echo running plink',
             'bcftools': 'echo running bcftools',
             'R': 'echo running R',
             'vcftools': 'echo running vcftools'}


subworkflow part1:
    snakefile: "part1.smk"

localrules: all

if config['check_vcf']:
    rule all:
        input:
            expand("log/{sample}_chr{chrom}.check.log",
                   sample=SAMPLE, chrom=CHROM)
else:
    rule all:
        input:
            expand("final/{sample}_chr{chrom}.vcf.gz",
                   sample=SAMPLE, chrom=CHROM)

if config['chr_callrate']:
    rule check_chr_callrate:
        input:
            "data/{sample}_chr{chrom}_preCallcheck.vcf.gz"
        output:
            "data/callrate/{sample}/chr{chrom}.sample_missingness.imiss"
        params:
            out = "data/callrate/{sample}/chr{chrom}.sample_missingness"
        shell:
            "{loads[vcftools]}; "
            "{com[vcftools]} --missing-indv --gzvcf {input} --out {params.out}"

    rule process_chr_callrate:
        input:
            expand("data/callrate/{{sample}}/chr{chrom}." +
                   "sample_missingness.imiss", chrom=CHROM)
        output:
            "data/callrate/{sample}/chrall.irem"
        params:
            dir = "data/callrate/{sample}"
        shell:
            "{loads[R]}; "
            "{com[R]} scripts/process_imiss.R {params.dir} 0.2"


def chunkfiles(wildcards):
    chrom = [str(x) for x in CHROM]
    fname = "data/callrate/{}.chunks.json".format(wildcards.sample)
    with open(fname, 'r') as chunkfile:
        chunks = json.load(chunkfile)
    chunks_proc = zip(chunks['chrom'], chunks['from'], chunks['through'])
    files = ['{}/chr{}_from{}_through{}'.format(wildcards.sample, ch, fr, to)
             for ch, fr, to in chunks_proc if ch in chrom]
    outname = "data/callrate/{}.sample_missingness.imiss"
    files = [outname.format(x) for x in files]
    return files


if config['chunk_callrate']:
    rule check_chunk_callrate:
        input:
            vcf = "data/{sample}_chrall_preCallcheck.vcf.gz"
        output:
            "data/callrate/{sample}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness.imiss"
        params:
            out = "data/callrate/{sample}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness",
            ranges = "--chr {chrom} --from-bp {range_from} --to-bp {range_through}"
        shell:
            "{loads[vcftools]}; "
            "{com[vcftools]} --missing-indv --gzvcf {input.vcf} {params.ranges} --out {params.out}"

    rule process_chunk_callrate:
        input: chunkfiles
        output:
            "data/callrate/{sample}/chrall.irem"
        params:
            dir = "data/callrate/{sample}"
        shell:
            "{loads[R]}; "
            "{com[R]} scripts/process_chunk_imiss.R {params.dir}"

rule apply_irem:
    input:
        irem = "data/callrate/{sample}/chrall.irem",
        vcf = "data/{sample}_chr{chrom}_preCallcheck.vcf.gz"
    output:
        "final/{sample}_chr{chrom}.vcf.gz"
    shell:
        "{loads[bcftools]}; "
        "{com[bcftools]} view -S ^{input.irem} -Oz -o {output} {input.vcf}; "
        "bcftools index -t {output}"

if config['check_vcf']:
    rule check_vcf:
        input:
            fasta = config["ref"],
            vcf = "final/{sample}_chr{chrom}.vcf.gz"
        output:
            "log/{sample}_chr{chrom}.check.log"
        params:
            py_env = config["checkEnv"]
        log: "log/{sample}_chr{chrom}"
        shell:
            "source activate {params.py_env}; "
            "python scripts/checkVCF.py -r {input.fasta} -o {log} {input.vcf}"
