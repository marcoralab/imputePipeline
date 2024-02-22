'''Snakefile for MIS preparation
   Version 1.0'''

BPLINK = ['bed', 'bim', 'fam']

#localrules: all, var_qc, subj_qc, split_to_vcf

# Pre-split QC

if 'maf' in config['preqc'] and config['preqc']['maf']:
    maf_cmd = '--maf ' + config['preqc']['maf']
else:
    maf_cmd = ''

rule var_qc:
    input: lambda wc: multiext(infiles[wc['cohort']]['fstem'], ".bim", ".bed", ".fam")
    output: temp(expand('{{outdir}}/plink/{{cohort}}_varqc.{ext}', ext=BPLINK))
    params:
        ins = lambda wc: infiles[wc['cohort']]['fstem'],
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
        keep_remove = sampfilt,
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

do_hwe = False
if 'hwe' in config['preqc'] and config['preqc']['hwe']:
    do_hwe = True
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
        plink = rules.hwe_qc.output if do_hwe else rules.subj_qc.output
    params:
        bim = rules.hwe_qc.params.out + '.bim' if do_hwe else rules.subj_qc.params.out + '.bim',
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
        mapping = '{outdir}/rename_chrom/{cohort}_mapping.txt'
    threads: 1
    resources:
        mem_mb = 256,
        time_min = 10
    container: 'docker://befh/flippyr:0.5.3'
    script: '../scripts/rename_chrom.py'

def get_chrname(wc):
    if wc['chrom'] in [str(x) for x in range(1, 23)]:
        return wc['chrom']
    if wc['chrom'] == 'X':
        return 'X,XY'
    if wc['chrom'] == 'Y':
        return 'Y'
    if wc['chrom'] == 'M':
        return 'MT'
    raise Exception('Bad chromosome')

rule split_to_vcf:  # Split plink files into chromosomes.
    input: rules.flippyr.output.plink
    params:
        ins = '{outdir}/plink/{cohort}_refmatched',
        out = '{outdir}/{cohort}.chr{chrom}_unsorted',
        chr = get_chrname
    output: temp('{outdir}/{cohort}.chr{chrom}_unsorted.vcf.gz')
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell: '''
plink --bfile {params.ins} --chr {params.chr} \
  --memory 16000 --real-ref-alleles \
  --recode vcf bgz --out {params.out} --silent
'''

rule rename_vcf_fromplink:
    input:
        vcf = rules.split_to_vcf.output,
        rename = rules.rename_chrom.output.mapping
    output: temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|M}_preCallcheck.vcf.gz')
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    shell:
        '''
bcftools annotate --rename-chrs -Oz -o {output} {input.rename} 
'''

def get_chrom(wc):
    if build == "b38":
        return f'chr{wc["chrom"]}'
    else:
        chrom = get_chrname(wc)
        return 'Y' if chrom == 'X,XY' else chrom

rule invcf_split:
    input:
        vcf = lambda wc: infiles[wc['cohort']]['file'],
        sampfilt = select_sampfilt_files(config)
    output: temp('{outdir}/prep/{cohort}_split_chr{chrom,[0-9XY]+|M}_varqc.vcf.gz')
    params:
        chrom = get_chrom,
        sf = sampfilt
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'    
    shell:
        '''
bcftools view{sf} -Oz --threads 2 --targets-file {params.chrom} -o {output}
'''

rule invcf_sampfilt:
    input:
        vcf = lambda wc: infiles[wc['cohort']]['file'],
        sampfilt = select_sampfilt_files(config)
    output: temp('{outdir}/prep/{cohort}_split_chr{chrom,[0-9XY]+|M}_varqc.vcf.gz')
    params:
        sf = sampfilt
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'    
    shell:
        '''
bcftools view{sf} -Oz --threads 2 -o {output}
'''

def input_filter_vars_vcf(wc):
    if infiles[wc['cohort']]['type'] == 'vcf':
        return rules.invcf_split.output
    elif sampfilt:
        return rules.invcf_sampfilt.output
    else:
        return infiles[wc['cohort']]['file']



rule filter_vars_vcf:
    input: input_filter_vars_vcf
    output: temp('{outdir}/prep/{cohort}_chr{chrom,[0-9XY]+|M}_varqc.vcf.gz')
    params:
        geno = config['preqc']['geno'],
        maf = f"MAF[0]>={config['preqc']['maf']} && " if 'maf' in config['preqc'] and config['preqc']['maf'] else ""
    threads: 3
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    shell:
        '''
bcftools +fill-tags {input} -- -t MAF,HWE,F_MISSING | \\
  bcftools view -i '{params.maf}F_MISSING<{params.geno}' -Oz --threads 2 -o {output}
'''

rule get_sampmiss_vcf:
    input: rules.filter_vars_vcf.output
    output: temp('{outdir}/prep/{cohort}_chr{chrom,[0-9XY]+|M}_sampmiss.smiss')
    params:
        out = '{outdir}/prep/{cohort}_chr{chrom,[0-9XY]+|M}_sampmiss'
    threads: 3
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell:
        '''
plink2 --vcf {input} --missing sample-only --out {params.out}
'''

rule proc_sampmiss_vcf:
    input: expand('{{outdir}}/prep/{{cohort}}_chr{chrom}_sampmiss.smiss', chrom = CHROM)
    output:
        smiss = '{outdir}/prep/{cohort}_allchr_sampmiss.smiss',
        ikeep = '{outdir}/prep/{cohort}_sampmiss.ikeep'
    params:
        mind = config['preqc']['mind']
    threads: 1
    resources:
        mem_mb = 4000,
        time_min = 30
    conda: '../envs/r.yaml'
    script: '../scripts/rule_proc_sampmiss_vcf.R'

rule apply_sampmiss_vcf:
    input:
        vcf = rules.filter_vars_vcf.output,
        ikeep = rules.proc_sampmiss_vcf.output.ikeep
    output: temp('{outdir}/prep/{cohort}_split_chr{chrom,[0-9XY]+|M}_sampmiss.vcf.gz')
    threads: 2
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'    
    shell:
        '''
bcftools view --samples-file {input.ikeep} -Oz --threads 2 -o {output}
'''

if do_hwe:
    rule filter_hwe_vcf:
        input: rules.apply_sampmiss_vcf.output
        output: temp('{outdir}/prep/{cohort}_chr{chrom,[0-9XY]+|M}_varqc.vcf.gz')
        params:
            hwe = config['preqc']['hwe']
        threads: 3
        resources:
            mem_mb = 8192,
            time_min = 30
        conda: '../envs/bcftools.yaml'
        shell:
            '''
bcftools view -i 'HWE>{params.hwe}' -Oz --threads 2 -o {output} {input}
'''

def input_sort_vcf_precallrate(wc):
    if infiles[wc['cohort']]['type'] == 'plink':
        return rules.rename_vcf_fromplink.output
    elif do_hwe:
        return rules.filter_hwe_vcf.output
    else:
        return rules.apply_sampmiss_vcf.output

rule sort_vcf_precallrate:
    input:
        vcf = rules.rename_vcf_fromplink.output
    output:
        vcf = temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|M}_preCallcheck.vcf.gz'),
        tbi = temp('{outdir}/{cohort}_chr{chrom,[0-9XY]+|M}_preCallcheck.vcf.gz.tbi')
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    shell:
        '''
bcftools sort {input} -Oz -o {output.vcf}
bcftools index -t {output.vcf}
'''