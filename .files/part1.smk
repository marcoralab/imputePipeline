'''Snakefile for MIS preparation
   Version 0.9'''

from scripts.parse_config import parser

configfile: "config.yaml"
shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]

CHROM, SAMPLE, INPATH, keep_command = parser(config)

localrules: all, var_qc, subj_qc, split_to_vcf

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

if isMinerva:
    com = {'flippyr': 'flippyr', 'plink': 'plink --keep-allele-order',
           'plink2': 'plink', 'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R',
           'king': 'king'}
    loads = {'flippyr': 'module load plink/1.90b6.7', 'plink': 'module load plink/1.90b6.7',
             'bcftools': 'module load bcftools/1.9',
             'king': 'module unload gcc; module load king/2.1.6',
             'R': ('module load R/3.5.3 pandoc/2.1.3 udunits/2.2.26; ',
                   'RSTUDIO_PANDOC=$(which pandoc)')}
else:
    com = {'flippyr': 'flippyr',
           'plink': 'plink --keep-allele-order', 'plink2': 'plink',
           'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R', 'king': 'king'}
    loads = {'flippyr': 'echo running flippyr', 'plink': 'echo running plink',
             'bcftools': 'echo running bcftools',  'R': 'echo running R',
             'king': 'echo running KING'}

if isMinerva:
    shell.prefix("PATH=" + config["anaconda"] + ":$PATH; ")

if config['chr_callrate']:
    rule all:
        input:
            expand("data/{sample}_chr{chrom}_preCallcheck.vcf.gz",
                   sample=SAMPLE, chrom=CHROM)
elif config['chunk_callrate']:
    rule all:
        input:
            expand("data/{sample}_chr{chrom}_preCallcheck.vcf.gz",
                   sample=SAMPLE, chrom=CHROM),
            expand("data/callrate/{sample}.chunks.json", sample=SAMPLE)
else:
    rule all:
        input:
            expand("final/{sample}_chr{chrom}.vcf.gz",
                   sample=SAMPLE, chrom=CHROM)

# Pre-split QC

rule var_qc:
    input:
        expand(INPATH + "{{sample}}.{ext}", ext=BPLINK)
    output:
        plink = expand("data/plink/{{sample}}_varqc.{ext}", ext=BPLINK)
    params:
        ins = INPATH + "{sample}",
        out = "data/plink/{sample}_varqc",
        geno = config["qc"]["geno"],
        hwe = config["qc"]["hwe"],
        maf = config["qc"]["maf"]
    threads: 1
    shell:
        "{loads[plink]}; "
        "{com[plink]} -bfile {params.ins} --geno {params.geno} --memory 128 "
        "--keep-allele-order --hwe {params.hwe} --maf {params.maf} "
        "--make-bed --out {params.out} --silent"

rule subj_qc:
    input:
        plink = rules.var_qc.output.plink
    params:
        ins = rules.var_qc.params.out,
        out = "data/plink/{sample}_indivqc",
        mind = config["qc"]["mind"],
        keep = keep_command
    output:
        expand("data/plink/{{sample}}_indivqc.{ext}", ext=BPLINK)
    threads: 1
    shell:
        "{loads[plink]}; "
        "plink -bfile {params.ins} --memory 128 "
        "--mind {params.mind} --remove samp.irem {params.keep}"
        "--make-bed --out {params.out} --silent"

rule flippyr:
    input:
        fasta = config["ref"],
        plink = rules.subj_qc.output
    params:
        out = "data/plink/{sample}",
        suff = "_refmatched"
    output:
        expand("data/plink/{{sample}}.{ext}",
               ext=["allele", "flip", "delete", "log", "log.tab"]),
        command = "data/plink/{sample}.runPlink",
        plink = expand("data/plink/{{sample}}_refmatched.{ext}", ext=BPLINK)
    run:
        import flippyr
        flippyr.writeFiles(input["fasta"], input["plink"][1], params["out"],
                           silent=False, plink=False, p_suff=params["suff"])
        shell("{plink}; bash {}".format(output["command"], plink=loads[plink]))

# Split, sort and compress

rule split_to_vcf:  # Split plink files into chromosomes.
    input:
        rules.flippyr.output.plink
    params:
        ins = "data/plink/{sample}_refmatched",
        out = "data/{sample}.chr{chrom}_unsorted",
        c = "{chrom}"
    output:
        "data/{sample}.chr{chrom}_unsorted.vcf"
    shell:
        "{loads[plink]}; "
        "{com[plink2]} -bfile {params.ins} --chr {params.c} "
        "--memory 256 --real-ref-alleles "
        "--recode vcf --out {params.out}"

rule all_to_vcf:
    input:
        rules.flippyr.output.plink
    params:
        ins = "data/plink/{sample}_refmatched",
        out = "data/{sample}_chrall_unsorted",
    output:
        "data/{sample}_chrall_unsorted.vcf"
    shell:
        "{loads[plink]}; "
        "plink -bfile {params.ins} "
        "--memory 256 --real-ref-alleles "
        "--recode vcf --out {params.out}"

rule compress_vcf_allchr:
    input:
        rules.all_to_vcf.output
    output:
        "data/{sample}_chrall_preCallcheck.vcf.gz"
    threads: 8
    shell:
        "{loads[bcftools]}; "
        "{com[bcftools]} sort {input} -Oz -o {output}; "
        "{com[bcftools]} index -t {output}"

rule make_chunk_yaml:
    input:
        rules.compress_vcf_allchr.output
    output:
        "data/callrate/{sample}.chunks.json"
    shell:
        "python scripts/fullchunker.py {input} {output}"

rule compress_vcf_precallrate:
    input:
        rules.split_to_vcf.output
    output:
        "data/{sample}_chr{chrom}_preCallcheck.vcf.gz"
    threads: 8
    shell:
        "{loads[bcftools]}; "
        "{com[bcftools]} sort {input} -Oz -o {output}; "
        "{com[bcftools]} index -t {output}"

rule compress_vcf_nocallrate:
    input:
        rules.split_to_vcf.output
    output:
        "final/{sample}_chr{chrom}.vcf.gz"
    threads: 8
    shell:
        "{loads[bcftools]}; "
        "{com[bcftools]} sort {input} -Oz -o {output}; "
        "{com[bcftools]} index -t {output}"
