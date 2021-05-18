'''Snakefile for MIS preparation
   Version 1.0'''

from scripts.parse_config import parser

configfile: 'config/config.yaml'
shell.executable('/bin/bash')

BPLINK = ['bed', 'bim', 'fam']

CHROM, SAMPLE, INPATH, keep_command = parser(config)

localrules: all, var_qc, subj_qc, split_to_vcf

# Pre-split QC

rule var_qc:
    input: expand(INPATH + '{{sample}}.{ext}', ext=BPLINK)
    output: expand('data/plink/{{sample}}_varqc.{ext}', ext=BPLINK)
    params:
        ins = INPATH + '{sample}',
        out = 'data/plink/{sample}_varqc',
        geno = config['qc']['geno'],
        hwe = config['qc']['hwe'],
        maf = config['qc']['maf']
    threads: 1
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order --bfile {params.ins} --memory 2560 \
  --geno {params.geno} --hwe {params.hwe} --maf {params.maf} \
  --make-bed --out {params.out} --silent
'''

rule subj_qc:
    input: rules.var_qc.output
    params:
        ins = rules.var_qc.params.out,
        out = 'data/plink/{sample}_indivqc',
        mind = config['qc']['mind'],
        keep = keep_command
    output: expand('data/plink/{{sample}}_indivqc.{ext}', ext=BPLINK)
    threads: 1
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order --bfile {params.ins} --memory 1280 \
  --mind {params.mind} --remove samp.irem {params.keep} \
  --make-bed --out {params.out} --silent
'''

rule flippyr:
    input:
        fasta = config['ref'],
        plink = rules.subj_qc.output
    params:
        out = 'data/plink/{sample}',
        suff = '_refmatched'
    output:
        expand('data/plink/{{sample}}.{ext}',
               ext=['allele', 'flip', 'delete', 'log', 'log.tab']),
        command = 'data/plink/{sample}.runPlink',
        plink = expand('data/plink/{{sample}}_refmatched.{ext}', ext=BPLINK)
    conda: 'envs/flippyr.yaml'
    shell: 'flippyr -o {params.out} --outputSuffix {params.suff} --plink'

# Split, sort and compress

rule split_to_vcf:  # Split plink files into chromosomes.
    input: rules.flippyr.output.plink
    params:
        ins = 'data/plink/{sample}_refmatched',
        out = 'data/{sample}.chr{chrom}_unsorted',
        c = '{chrom}'
    output: temp('data/{sample}.chr{chrom}_unsorted.vcf.bgz')
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --bfile {params.ins} --chr {params.c} --memory 256 --real-ref-alleles \
  --recode vcf bgz --out {params.out}
'''

rule sort_vcf_precallrate:
    input: rules.split_to_vcf.output
    output:
        vcf = temp('data/{sample}_chr{chrom}_preCallcheck.vcf.gz'),
        tbi = temp('data/{sample}_chr{chrom}_preCallcheck.vcf.gz.tbi')
    threads: 8
    conda: 'envs/bcftools.yaml'
    shell:
        '''
bcftools sort {input} -Oz -o {output.vcf}
bcftools index -t {output.vcf}
'''
