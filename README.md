# imputePrep

A Snakemake pipeline to prepare binary PLINK files for the Michigan Imputation Server.

## Requirements:

 * Anaconda or Miniconda python
   * (Alternatively make your own Python2 environment and copy or symlink to the `scripts/checkEnv` folder of your install location. You will also have to manually copy the other files.)
 * Python3 (The future is now!)
 * Snakemake (Please install first.)

## Usage:

### Config files:

There are two config files with settings for imputation prep:

 * *config.yaml* contains settings for the QC and preparation.
 * *cluster.yaml* contains settings for cluster execution.

Please review those files before starting.

### Data-file location:

By default, the pipeline looks for PLINK filesets in the base directory of your installation. You can choose a different destination by editing `config.yaml`.
