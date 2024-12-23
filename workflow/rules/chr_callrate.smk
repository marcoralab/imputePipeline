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
    conda: '../envs/bcftools.yaml'
    shell: '''
vcf={input.vcf}
bcftools stats -S <(bcftools query -l $vcf) $vcf | \
  awk 'BEGIN {{FS=OFS="\t"; print "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS"}} \
       $1 == "SN" && $3 == "number of records:" {{n=$4}} \
       $1 == "PSC" {{print $3,n,0,$14,$14/n}}' > {output}
'''

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
    conda: '../envs/r.yaml'
    script: '../scripts/process_imiss.R {params.dir}'
