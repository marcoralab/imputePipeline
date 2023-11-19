'''Snakefile for MIS preparation
   Version 1.0'''

BPLINK = ['bed', 'bim', 'fam']

#localrules: all, var_qc, subj_qc, split_to_vcf

# Pre-split QC

if config['preqc']['maf']:
    maf_cmd = '--maf ' + config['preqc']['maf']
else:
    maf_cmd = ''

rule var_qc:
    input: multiext(INPATH + '/{cohort}', ".bim", ".bed", ".fam")
    output: temp(expand('{{outdir}}/plink/{{cohort}}_varqc.{ext}', ext=BPLINK))
    params:
        ins = INPATH + '/{cohort}',
        out = '{outdir}/plink/{cohort}_varqc',
        geno = config['preqc']['geno'],
        maf = maf_cmd
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
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
        mind = config['preqc']['mind'],
        keep_remove = keep_remove_command,
    output: temp(expand('{{outdir}}/plink/{{cohort}}_indivqc.{ext}', ext=BPLINK))
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell:
        '''
plink --keep-allele-order \
  --bfile {params.ins} --memory 8192 \
  --mind {params.mind} {params.keep_remove} \
  --make-bed --out {params.out} --silent
'''

if config['preqc']['hwe']:
    rule hwe_qc:
        input: rules.subj_qc.output
        output: temp(expand('{{outdir}}/plink/{{cohort}}_hwe.{ext}', ext=BPLINK))
        params:
            ins = rules.subj_qc.params.out,
            out = '{outdir}/plink/{cohort}_hwe',
            hwe = config['preqc']['hwe'],
        threads: 2
        resources:
            mem_mb = 8192,
            time_min = 30
        conda: '../envs/plink.yaml'
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
        plink = rules.hwe_qc.output if config['preqc']['hwe'] else rules.subj_qc.output
    params:
        bim = rules.hwe_qc.params.out + '.bim' if config['preqc']['hwe'] else rules.subj_qc.params.out + '.bim',
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

rule rename_chrom:
    input:
        fasta = config['ref'],
        bim = '{outdir}/plink/{cohort}_refmatched.bim'
    output:
        json = '{outdir}/rename_chrom/{cohort}_mapping.json',
        bim = '{outdir}/rename_chrom/{cohort}_refmatched_renamed.bim'
    threads: 1
    resources:
        mem_mb = 256,
        time_min = 10
    container: 'docker://befh/flippyr:0.5.3'
    script: '../scripts/rename_chrom.py'

rule get_chrname:  # Split plink files into chromosomes.
    input: rules.rename_chrom.output.json
    output: temp('{outdir}/{cohort}_{chrom}.name')
    threads: 1
    resources:
        mem_mb = 128,
        time_min = 1
    run:
        import json

        with open(input[0], 'r') as f:
            m = json.load(f)

        if int(m['build']) == 37 and wildcards.chrom == 'M':
            chr_ = 'MT'
        elif int(m['build']) == 37:
            chr_ = wildcards.chrom
        else:
            chr_ = m['map'][wildcards.chrom]

        with open(output[0], 'w') as f:
            print(chr_, file=f)

rule sort_plink:  # Split plink files into chromosomes.
    input:
        fileset = rules.flippyr.output.plink,
        bim = rules.rename_chrom.output.bim
    params:
        ins = '{outdir}/plink/{cohort}_refmatched',
        out = '{outdir}/rename_chrom/{cohort}_flip-rename-sort'
    output:
        temp(multiext('{outdir}/rename_chrom/{cohort}_flip-rename-sort',
                      '.bed', '.bim', '.fam'))
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell: '''
plink --bfile {params.ins} --bim {input.bim} --real-ref-alleles \
  --make-bed --out {params.out}
'''

rule split_to_vcf:  # Split plink files into chromosomes.
    input:
        fileset = rules.sort_plink.output,
        chrname = rules.get_chrname.output
    params:
        ins = '{outdir}/rename_chrom/{cohort}_flip-rename-sort',
        out = '{outdir}/{cohort}.chr{chrom}_unsorted'
    output: temp('{outdir}/{cohort}.chr{chrom}_unsorted.vcf.gz')
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell: '''
plink --bfile {params.ins} --chr "$(cat {input.chrname})" \
  --memory 16000 --real-ref-alleles \
  --recode vcf bgz --out {params.out} --silent
'''

rule sort_vcf_precallrate:
    input: rules.split_to_vcf.output
    output:
        vcf = temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|MT}_preCallcheck.vcf.gz'),
        tbi = temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|MT}_preCallcheck.vcf.gz.tbi')
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    shell:
        '''
bcftools sort -Oz -o {output.vcf} {input}
bcftools index -t {output.vcf}
'''
