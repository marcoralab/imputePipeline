suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

inputs <- snakemake@input

read_smiss <- function(fname) {
  read_tsv(
    fname,
    col_types = cols(`#IID` = "c",
                     MISSING_CT = "i",
                     OBS_CT = "i",
                     .default = "-")) |>
    rename(ID = "#IID")
}

add_smiss <- function(smiss, fname) {
  if (is.null(smiss)) {
    return(read_smiss(fname))
  }
  smiss |>
    full_join(read_smiss(fname), by = "ID") |>
    mutate(across(where(is.numeric), ~replace_na(., 0)),
           MISSING_CT = MISSING_CT.x + MISSING_CT.y,
           OBS_CT = OBS_CT.x + OBS_CT.y) |>
    select(ID, MISSING_CT, OBS_CT)
}

smiss <- NULL
for (input in inputs) {
  smiss <- add_smiss(smiss, input)
}

smiss_tot <- smiss |>
  mutate(F_MISS = MISSING_CT / OBS_CT) |>
  write_tsv(snakemake@output[['smiss']])
  
smiss_tot |>
  filter(F_MISS < as.double(snakemake@params[['mind']])) |>
  pull(ID) |>
  write_tsv(snakemake@output[['ikeep']])