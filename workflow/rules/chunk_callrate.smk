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

rule all_to_vcf:
    input: rules.flippyr.output.plink
    params:
        ins = '{outdir}/plink/{cohort}_refmatched',
        out = '{outdir}/{cohort}_chrall_unsorted',
    output: temp('{outdir}/{cohort}_chrall_unsorted.vcf.gz')
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/plink.yaml'
    shell:
        '''
plink --bfile {params.ins} --memory 4096 --real-ref-alleles \
  --recode vcf bgz --out {params.out}
'''

rule sort_vcf_allchr:
    input: rules.all_to_vcf.output
    output:
        vcf = temp('{outdir}/{cohort}_chrall_preCallcheck.vcf.gz'),
        tbi = temp('{outdir}/{cohort}_chrall_preCallcheck.vcf.gz.tbi')
    params:
        tempdir = "{outdir}/temp/{cohort}"
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 120
    conda: '../envs/bcftools.yaml'
    shell:
        '''
mkdir -p {params.tempdir}
bcftools sort -Oz -o {output.vcf} \
  --max-mem 64000M -T {params.tempdir} {input}
bcftools index -t {output.vcf}
'''

checkpoint make_chunk_yaml:
    input:
        vcf = rules.sort_vcf_allchr.output.vcf,
        tbi = rules.sort_vcf_allchr.output.tbi
    output: '{outdir}/callrate/{cohort}.chunks.json'
    threads: 4
    resources:
        mem_mb = 8192,
        time_min = 30
    conda: '../envs/chunking.yaml'
    script: 'scripts/fullchunker.py'

rule check_chunk_callrate:
    input:
        vcf = rules.sort_vcf_allchr.output.vcf,
        tbi = rules.sort_vcf_allchr.output.tbi
    output:
        '{outdir}/callrate/{cohort}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness.imiss'
    params:
        out = '{outdir}/callrate/{cohort}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness',
        ranges = '--chr {chrom} --from-bp {range_from} --to-bp {range_through}'
    threads: 1
    resources:
        mem_mb = 5200,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    shell:
        '''
vcftools --missing-indv --gzvcf {input.vcf} {params.ranges} --out {params.out}
'''
# import ipdb; ipdb.set_trace()
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
        time_min = 30
    conda: '../envs/r.yaml'
    script: 'scripts/process_chunk_imiss.R'
