rule check_chr_callrate:
    input:
        vcf = '{outdir}/{cohort}_chr{chrom}_preCallcheck.vcf.gz',
        tbi = '{outdir}/{cohort}_chr{chrom}_preCallcheck.vcf.gz.tbi'
    output: '{outdir}/callrate/{cohort}/chr{chrom}.sample_missingness.imiss'
    params:
        out = '{outdir}/callrate/{cohort}/chr{chrom}.sample_missingness'
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: 'envs/bcftools.yaml'
    shell: 'vcftools --missing-indv --gzvcf {input.vcf} --out {params.out}'

rule process_chr_callrate:
    input:
        expand('{{outdir}}/callrate/{{cohort}}/chr{chrom}.' +
               'sample_missingness.imiss', chrom=CHROM)
    output: '{outdir}/callrate/{cohort}/chrall.irem'
    params:
        dir = '{outdir}/callrate/{cohort}',
        threshold = 0.2
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: 'envs/r.yaml'
    script: 'scripts/process_imiss.R {params.dir}'
