def chunkfiles(wildcards):
    chrom = [str(x) for x in CHROM]
    fname = checkpoints.make_chunk_yaml.get(sample=wildcards.sample).output[0]
    with fname.open() as chunkfile:
        chunks = json.load(chunkfile)
    chunks_proc = zip(chunks['chrom'], chunks['from'], chunks['through'])
    files = ['{}/chr{}_from{}_through{}'.format(wildcards.sample, ch, fr, to)
             for ch, fr, to in chunks_proc if ch in chrom]
    outname = 'data/callrate/{}.sample_missingness.imiss'
    files = [outname.format(x) for x in files]
    return files

rule all_to_vcf:
    input: rules.flippyr.output.plink
    params:
        ins = 'data/plink/{sample}_refmatched',
        out = 'data/{sample}_chrall_unsorted',
    output: temp('data/{sample}_chrall_unsorted.vcf.gz')
    conda: 'envs/plink.yaml'
    shell:
        '''
plink --bfile {params.ins} --memory 256 --real-ref-alleles \
  --recode vcf bgz --out {params.out}
'''

rule sort_vcf_allchr:
    input: rules.all_to_vcf.output
    output:
        vcf = temp('data/{sample}_chrall_preCallcheck.vcf.gz'),
        tbi = temp('data/{sample}_chrall_preCallcheck.vcf.gz.tbi')
    threads: 8
    conda: 'envs/bcftools.yaml'
    shell:
        '''
bcftools sort {input} -Oz -o {output.vcf}
bcftools index -t {output.vcf}
'''

checkpoint make_chunk_yaml:
    input:
        vcf = rules.sort_vcf_allchr.output.vcf,
        tbi = rules.sort_vcf_allchr.output.tbi
    output: 'data/callrate/{sample}.chunks.json'
    script: 'scripts/fullchunker.py'

rule check_chunk_callrate:
    input:
        vcf = rules.sort_vcf_allchr.output.vcf,
        tbi = rules.sort_vcf_allchr.output.tbi
    output:
        'data/callrate/{sample}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness.imiss'
    params:
        out = 'data/callrate/{sample}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness',
        ranges = '--chr {chrom} --from-bp {range_from} --to-bp {range_through}'
    conda: 'envs/bcftools.yaml'
    shell:
        '''
vcftools --missing-indv --gzvcf {input.vcf} {params.ranges} --out {params.out}'
'''

rule process_chunk_callrate:
    input: chunkfiles
    output: 'data/callrate/{sample}/chrall.irem'
    params:
        dir = 'data/callrate/{sample}'
    conda: 'envs/r.yaml'
    shell: 'Rscript scripts/process_chunk_imiss.R {params.dir}'
