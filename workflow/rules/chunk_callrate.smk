def chunkfiles(wc):
    chrom = [str(x) for x in CHROM]
    fname = checkpoints.make_chunk_yaml.get(
        cohort=wc.cohort, outdir=wc.outdir).output[0]
    with fname.open() as chunkfile:
        chunks = json.load(chunkfile)
    chunks_proc = zip(chunks['chrom'], chunks['from'], chunks['through'])
    ranges = [f'chr{ch}_from{fr}_through{to}'
              for ch, fr, to in chunks_proc if ch in chrom]
    files = [
        f'{wc.outdir}/callrate/{wc.cohort}/{range}.sample_missingness.imiss'
        for range in ranges]
    return files

### This is slow and inefficient. Better to make chunks per chromosome then merge the json files

rule cat_chroms:
    input:
        vcf = expand('{{outdir}}/{{cohort}}_chr{chrom}_preCallcheck.vcf.gz', chrom=CHROM),
        tbi = expand('{{outdir}}/{{cohort}}_chr{chrom}_preCallcheck.vcf.gz.tbi', chrom=CHROM)
    output:
        vcf = '{outdir}/{cohort}_allchroms_preCallcheck.vcf.gz',
        tbi = '{outdir}/{cohort}_allchroms_preCallcheck.vcf.gz.tbi'
    threads: 1
    resources:
        mem_mb = 5200,
        time_min = 600
    conda: '../envs/bcftools.yaml'
    shell: '''
bcftools concat {input.vcf} -Oz -o {output.vcf}
bcftools index -tf {output.vcf}
'''

checkpoint make_chunk_yaml:
    input:
        vcf = rules.cat_chroms.output.vcf,
        tbi = rules.cat_chroms.output.tbi
    output: '{outdir}/callrate/{cohort}.chunks.json'
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 60
    conda: '../envs/chunking.yaml'
    script: '../scripts/fullchunker.py'

rule check_chunk_callrate:
    input:
        vcf = rules.sort_vcf_precallrate.output.vcf,
        tbi = rules.sort_vcf_precallrate.output.tbi
    output:
        '{outdir}/callrate/{cohort}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness.imiss'
    params:
        out = '{outdir}/callrate/{cohort}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness',
        ranges = '--chr {chrom} --from-bp {range_from} --to-bp {range_through}'
    threads: 1
    resources:
        mem_mb = 5200,
        time_min = 600
    conda: '../envs/bcftools.yaml'
    shell:
        '''
vcftools --missing-indv --gzvcf {input.vcf} {params.ranges} --out {params.out}
'''

rule process_chunk_callrate:
    input: chunkfiles
    output: '{outdir}/callrate/{cohort}/chrall.irem'
    params:
        dir = '{outdir}/callrate/{cohort}',
        threshold = 0.5,
        chunk_variant_count_min = 50
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 60
    conda: '../envs/r.yaml'
    script: '../scripts/process_chunk_imiss.R'
