# imputePrep

A Snakemake pipeline to prepare binary PLINK files for the Michigan Imputation Server.

## Requirements:

 * Anaconda or Miniconda python
 * Singularity
 * Python3 (The future is now!)
 * Snakemake

All other reqirements will be installed locally by the pipeline at runtime.

## Usage:

### Configuration:

`config/config.yaml` contains settings for QC and preparation, submission and download, and post-imputation processing.

Please review and edit the configuration file before starting.

#### Data-file location:

By default, the pipeline looks for PLINK filesets in the base directory of your installation. You can choose a different destination by editing `directory`. Cohort names are infered from the stem of the plink filesets.

If you do not want all filesets in the `directory` to be processed and imputed, then you can include `COHORT:` or `SAMPLES:` and a list of cohorts in the config file like this:

```yaml
COHORT:
  - group1
  - study2
```

or

```yaml
COHORT: [group1, study2]
```

#### Output directory:

The output files for the pipeline will be stored in the directory specified by `out_dir`. This directory *must* exist.

#### Output file choices

Output file choices can be specified under `outputs`, and you can delete or comment out any list item you do not want. The options follow:

  - `stat_report` A report of imputation quality
  - `vcf_bycohort` Bgzipped VCF files for each cohort
  - `vcf_merged` A bgzipped VCF with all cohorts
  - `bgen_bycohort` Binary oxford filesets for each cohort
  - `bgen_merged` A binary oxford fileset with all cohorts
  - `plink_bycohort` Binary plink filesets for each cohort
  - `plink_merged` A binary plink fileset with all cohorts

#### Imputation configuration

This pipeline can impute using the [NIH imputatation server](https://imputation.biodatacatalyst.nhlbi.nih.gov/) and the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!). You must have an account and API key with any server you plan to use.

Imputation settings are stored under `imputation:` as a hash/dictionary of cohort names, with a hash of settings underneath. You can provide default settings with the hash key `default`, and/or override those settings with your individual cohort names.

The hash of settings follows the API specifications [here](https://imputationserver.readthedocs.io/en/latest/api/#job-submission-for-whole-genome-imputation) and [here](https://topmedimpute.readthedocs.io/en/stable/api/#job-submission), in addition to your API key, stored under `token`. `files` and `password` are automatically set. Additionally, if you are planning to use HRC or TOPMed, and you just provide the refpanel or server along with your API key, best-practice defaults will be used.

Options for `server` are `NIH` or `Michigan`. Case does not matter.

For each cohort specified, the pipeline will override the defaults, with the specified settings. If those are unchanged from the default, they do not need to be edited. 

Here is an example:

```yaml
imputation:
  default:
    token: token_here
    refpanel: topmed-r2
  study2:
    token: other_token_here
    server: Michigan
    refpanel: gasp-v2
    population: ASN
```

This will run study2 on the Michigan server using the GAsP panel with the ASN population, and all other cohorts with TOPMed using all populations.

#### Chromosome selection

Select the chromosomes you want to upload by editing `chroms`. Separate ranges with ":" and listed chromosomes with ",". You can put both in the same string. Use M for mitochondra.

Options are 1:22,X,Y,M

#### Reference and genome build

The pipeline will reference-align the input files before imputation using the fasta file specified under `ref:`. This fasta MUST be in the same genome build as the input. Input genome builds must match the builds supported by the chosen server, and you must specify the build under `imputation` if it does not match the default for the server.

This pipeline supports all builds for the downloaded imputed files. Be aware that TOPMed imputation returns GRCh38 genotypes.

#### QC

The pipeline can both do QC in preparation for upload and filter the files provided by the servers on a per-cohort basis. Pre-imputation QC is under `preqc:` and post-imputation QC is under `postqc:`. The options are documented in the configuration file.

##### Sample fitering to remove low-callrate subjects

The imputation server imputes in 10kBase chunks, removing chunks that have below a 50% callrate in *any* subject. To avoid this, we can filter such that no chromosome has above 20% missingness (`chr_callrate: True`), or that no imputation chunk with at least 50 variants has a missingness greater than 50% by removing subjects who violate the criterion. You can skip both checks by setting both options to false, but you can only perform a maximum of one of these checks. The chunk check is recomended.

#### Sample inclusion or excusion by name

You can also include OR exclude samples by a named list file using one of the following options:

- `include_samp:` filename of samples to include (plink --keep)
- `exclude_samp:` filename of samples to exclude (plink --remove)

`include_samp` accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. `exclude_samp` does the same for all *listed* samples.

### Running

You must run the pipeline with Anaconda environments and Singularity enabled. You can either use a Snakemake profile to do so or run the pipeline with the following command:

```bash
snakemake --use-conda --use-singularity
```

Make sure your cluster supports the amount of memory required for report generation (528 GiB) and that the memory is being properly requested. If that is not the case, you can edit the resource requirements on the rule `stats` in the [postImpute](https://github.com/marcoralab/postImpute) pipeline and modify the modularization lines in `workflow/Snakefile` in this pipeline. We have included a `lsf.yaml` file that ensures those resources are available within this pipeline.


