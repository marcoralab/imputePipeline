#!/usr/bin/env Rscript

message("Loading Packages")

packload <- function(...) { # partially from pacman
  packages <- as.character(match.call(expand.dots = FALSE)[[2]])
  pl <- function(x) {
    suppressPackageStartupMessages(require(x, character.only = T))
  }
  no_output <- sapply(packages, pl)
}

packload(dplyr, readr, tibble, magrittr, stringr)

path <- path.expand(snakemake@params[["dir"]])
threshold <- as.double(snakemake@params[["threshold"]])
chunk_variant_count_min <- as.integer(
  snakemake@params[["chunk_variant_count_min"]])
f_list <- snakemake@input

#f_list <- list.files(path = path, pattern = "*.sample_missingness.imiss") %>%
#  paste0(path, "/", .)

message("Reading Missingness")

df <- tibble(INDV = character())
for (imiss in f_list) {
  chrom <- str_match(basename(imiss), "(chr.+)\\.sam")[1, 2]
  df2 <- read_tsv(imiss, col_types = cols(
    .default = "i",
    INDV = "c",
    F_MISS = "d"))
  n <- max(df2$N_DATA)
  if (n > chunk_variant_count_min) {
    df <- df2 %>%
      select(INDV, !!chrom := F_MISS) %>%
      left_join(df, by = "INDV")
  }
}

message("Processing Missingness")

df %<>%
  select_if(function(col) max(col) >= threshold) %>%
  mutate(maxs = apply(Filter(is.numeric, .), 1, max)) %>%
  filter(maxs >= threshold) %>%
  arrange(desc(maxs))

write_tsv(df, path = paste0(path, "/chrall.imiss"))
write_lines(df$INDV, paste0(path, "/chrall.irem"))
