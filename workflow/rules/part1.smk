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
    input:
        fileset = rules.flippyr.output.plink,
        chrname = rules.get_chrname.output
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

rule sort_vcf_precallrate:
    input:
        vcf = rules.split_to_vcf.output,
        rename = rules.rename_chrom.output.mapping
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
bcftools sort {input.vcf} | \
  bcftools annotate --rename-chrs {input.rename} -Oz -o {output.vcf} 
bcftools index -t {output.vcf}
'''
