'''Snakefile for MIS preparation
   Version 1.0'''

from scripts.parse_config import parser

configfile: 'config/config.yaml'
shell.executable('/bin/bash')

BPLINK = ['bed', 'bim', 'fam']

CHROM, SAMPLE, INPATH, keep_command = parser(config)

#localrules: all, var_qc, subj_qc, split_to_vcf

# Pre-split QC

if config['qc']['maf']:
    maf_cmd = '--maf ' + config['qc']['maf']
else:
    maf_cmd = ''

rule var_qc:
    input: expand(INPATH + '{{sample}}.{ext}', ext=BPLINK)
    output: temp(expand('data/plink/{{sample}}_varqc.{ext}', ext=BPLINK))
    params:
        ins = INPATH + '{sample}',
        out = 'data/plink/{sample}_varqc',
        geno = config['qc']['geno'],
        maf = maf_cmd
    threads: 1
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order --allow-extra-chr --chr 1-26,X,Y,XY,MT \
  --bfile {params.ins} --memory 8192 \
  --geno {params.geno} {params.maf} \
  --make-bed --out {params.out} --silent
'''

rule subj_qc:
    input: rules.var_qc.output
    params:
        ins = rules.var_qc.params.out,
        out = 'data/plink/{sample}_indivqc',
        mind = config['qc']['mind'],
        keep = keep_command
    output: temp(expand('data/plink/{{sample}}_indivqc.{ext}', ext=BPLINK))
    threads: 1
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order \
  --bfile {params.ins} --memory 8192 \
  --mind {params.mind} --remove samp.irem {params.keep} \
  --make-bed --out {params.out} --silent
'''

if config['qc']['hwe']:
    rule hwe_qc:
        input: rules.subj_qc.output
        output: temp(expand('data/plink/{{sample}}_hwe.{ext}', ext=BPLINK))
        params:
            ins = rules.subj_qc.params.out,
            out = 'data/plink/{sample}_hwe',
            hwe = config['qc']['hwe'],
        threads: 1
        conda: 'envs/plink.yaml'
        shell:
            '''
plink --keep-allele-order \
  --bfile {params.ins} --memory 8192 \
  --hwe {params.hwe} 'midp' \
  --make-bed --out {params.out} --silent
'''

rule flippyr:
    input:
        fasta = config['ref'],
        plink = rules.hwe_qc.output if config['qc']['hwe'] else rules.subj_qc.output
    params:
        bim = rules.hwe_qc.params.out + '.bim' if config['qc']['hwe'] else rules.subj_qc.params.out + '.bim',
        out = 'data/plink/{sample}',
        suff = '_refmatched'
    output:
        expand('data/plink/{{sample}}.{ext}',
               ext=['allele', 'flip', 'delete', 'log', 'log.tab']),
        command = 'data/plink/{sample}.runPlink',
        plink = temp(expand('data/plink/{{sample}}_refmatched.{ext}', ext=BPLINK))
    conda: 'envs/flippyr.yaml'
    shell: '''
flippyr -o {params.out} --outputSuffix {params.suff} --plink \
  {input.fasta} {params.bim}
'''

# Split, sort and compress

rule split_to_vcf:  # Split plink files into chromosomes.
    input: rules.flippyr.output.plink
    params:
        ins = 'data/plink/{sample}_refmatched',
        out = 'data/{sample}.chr{chrom}_unsorted',
        c = '{chrom}'
    output: temp('data/{sample}.chr{chrom}_unsorted.vcf.gz')
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --bfile {params.ins} --chr {params.c} --memory 8192 --real-ref-alleles \
  --recode vcf bgz --out {params.out} --silent
'''

rule sort_vcf_precallrate:
    input: rules.split_to_vcf.output
    output:
        vcf = temp('data/{sample}_chr{chrom,[0-9XY]+|MT}_preCallcheck.vcf.gz'),
        tbi = temp('data/{sample}_chr{chrom}_preCallcheck.vcf.gz.tbi')
    threads: 8
    conda: 'envs/bcftools.yaml'
    shell:
        '''
bcftools sort -Oz -o {output.vcf} {input}
bcftools index -t {output.vcf}
'''
