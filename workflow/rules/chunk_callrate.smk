def chunkfiles(wc):
    chrom = [str(x) for x in CHROM]
    fname = checkpoints.make_chunk_yaml.get(
        cohort=wc.cohort, outdir=wc.outdir).output[0]
    with fname.open() as chunkfile:
        chunks = json.load(chunkfile)
    chunks_proc = zip([re.sub("chr", "", x) for x in chunks['chrom']],
                      chunks['from'], chunks['through'])
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

#TODO rewrite to do per chromosome 

checkpoint make_chunk_yaml:
    input:
        vcf = rules.cat_chroms.output.vcf,
        tbi = rules.cat_chroms.output.tbi
    output: '{outdir}/callrate/{cohort}.chunks.json'
    threads: 4
    resources:
        mem_mb = 32768,
        time_min = 60
    conda: '../envs/chunking.yaml'
    script: '../scripts/fullchunker.py'

def check_chunk_callrate_set_ranges(wc):
    fname = checkpoints.rename_chrom.get(cohort=wc.cohort, outdir=wc.outdir).output.json
    with open(fname, "r") as f:
        build = json.load(f)['build']
    if build == 37:
        return f"{wc.chrom}:{wc.range_from}-{wc.range_through}"
    elif build == 38:
        return f"chr{wc.chrom}:{wc.range_from}-{wc.range_through}"
    else:
        raise ValueError(f"Invalid build {build}")


def check_chunk_callrate_buildfile(wc):
    fname = checkpoints.rename_chrom.get(cohort=wc.cohort, outdir=wc.outdir)
    return str(fname.output[0])

rule check_chunk_callrate:
    input:
        vcf = rules.sort_vcf_precallrate.output.vcf,
        tbi = rules.sort_vcf_precallrate.output.tbi,
        json = check_chunk_callrate_buildfile
    output:
        '{outdir}/callrate/{cohort}/chr{chrom}_from{range_from}_through{range_through}.sample_missingness.imiss'
    params:
        ranges = check_chunk_callrate_set_ranges
    threads: 1
    resources:
        mem_mb = 5200,
        time_min = 600
    conda: '../envs/bcftools.yaml'
    shell:
        '''
vcf={input.vcf} 
bcftools stats -S <(bcftools query -l $vcf) -r {params.ranges} $vcf | \
  awk 'BEGIN {{FS=OFS="\t"; print "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS"}} \
       $1 == "SN" && $3 == "number of records:" {{n=$4}} \
       $1 == "PSC" {{print $3,n,0,$14,$14/n}}' > {output}
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
        mem_mb = 32768,
        time_min = 60
    conda: '../envs/r.yaml'
    script: '../scripts/process_chunk_imiss.R'
