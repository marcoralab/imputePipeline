'''Snakefile for MIS preparation
   Version 1.0'''

from scripts.parse_config_imputePrep import parser

configfile: 'config/config.yaml'
shell.executable('/bin/bash')

BPLINK = ['bed', 'bim', 'fam']

CHROM, COHORT, INPATH, keep_command = parser(config)

#localrules: all, var_qc, subj_qc, split_to_vcf

# Pre-split QC

if config['qc']['maf']:
    maf_cmd = '--maf ' + config['qc']['maf']
else:
    maf_cmd = ''

rule var_qc:
    # input: expand(INPATH + '{{cohort}}.{ext}', ext=BPLINK)
    input: multiext('{cohort}', ".bim", ".bed", ".fam")
    output: temp(expand('{{outdir}}/plink/{{cohort}}_varqc.{ext}', ext=BPLINK))
    params:
        # ins = INPATH + '{cohort}',
        ins = '{cohort}',
        out = '{outdir}/plink/{cohort}_varqc',
        geno = config['qc']['geno'],
        maf = maf_cmd
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order --allow-extra-chr --chr 1-26,X,Y,XY,MT \
  --bfile {params.ins} --memory 8192 \
  --geno {params.geno} {params.maf} \
  --make-bed --out {params.out} --silent
'''

rule subj_qc:
    input: expand('{{outdir}}/plink/{{cohort}}_varqc.{ext}', ext=BPLINK)
    params:
        ins = rules.var_qc.params.out,
        out = '{outdir}/plink/{cohort}_indivqc',
        mind = config['qc']['mind'],
        keep = keep_command
    output: temp(expand('{{outdir}}/plink/{{cohort}}_indivqc.{ext}', ext=BPLINK))
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
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
        output: temp(expand('{{outdir}}/plink/{{cohort}}_hwe.{ext}', ext=BPLINK))
        params:
            ins = rules.subj_qc.params.out,
            out = '{outdir}/plink/{cohort}_hwe',
            hwe = config['qc']['hwe'],
        threads: 2
        resources:
            mem_mb = 8192,
            time_min = 30
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
        out = '{outdir}/plink/{cohort}',
        suff = '_refmatched'
    output:
        expand('{{outdir}}/plink/{{cohort}}.{ext}',
               ext=['allele', 'flip', 'delete', 'log', 'log.tab']),
        command = '{outdir}/plink/{cohort}.runPlink',
        plink = temp(expand('{{outdir}}/plink/{{cohort}}_refmatched.{ext}', ext=BPLINK))
    threads: 1
    resources:
        mem_mb = 8192,
        walltime = '48:00'
    container: 'docker://befh/flippyr:0.5.3'
    shell: '''
flippyr -o {params.out} --outputSuffix {params.suff} --plink \
  {input.fasta} {params.bim}
'''

# Split, sort and compress

rule split_to_vcf:  # Split plink files into chromosomes.
    input: rules.flippyr.output.plink
    params:
        ins = '{outdir}/plink/{cohort}_refmatched',
        out = '{outdir}/{cohort}.chr{chrom}_unsorted',
        c = '{chrom}'
    output: temp('{outdir}/{cohort}.chr{chrom}_unsorted.vcf.gz')
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --bfile {params.ins} --chr {params.c} --memory 16000 --real-ref-alleles \
  --recode vcf bgz --out {params.out} --silent
'''

rule sort_vcf_precallrate:
    input: rules.split_to_vcf.output
    output:
        vcf = temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|MT}_preCallcheck.vcf.gz'),
        tbi = temp('{outdir}/{cohort}_chr{chrom}_preCallcheck.vcf.gz.tbi')
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: 'envs/bcftools.yaml'
    shell:
        '''
bcftools sort --threads {threads} -Oz -o {output.vcf} {input}
bcftools index -t {output.vcf}
'''
