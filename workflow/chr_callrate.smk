rule check_chr_callrate:
    input:
        vcf = 'data/{sample}_chr{chrom}_preCallcheck.vcf.gz',
        tbi = 'data/{sample}_chr{chrom}_preCallcheck.vcf.gz.tbi'
    output: 'data/callrate/{sample}/chr{chrom}.sample_missingness.imiss'
    params:
        out = 'data/callrate/{sample}/chr{chrom}.sample_missingness'
    conda: 'envs/bcftools.yaml'
    shell: 'vcftools --missing-indv --gzvcf {input.vcf} --out {params.out}'

rule process_chr_callrate:
    input:
        expand('data/callrate/{{sample}}/chr{chrom}.' +
               'sample_missingness.imiss', chrom=CHROM)
    output: 'data/callrate/{sample}/chrall.irem'
    params:
        dir = 'data/callrate/{sample}',
        threshold = 0.2
    conda: 'envs/r.yaml'
    script: 'scripts/process_imiss.R {params.dir}'
