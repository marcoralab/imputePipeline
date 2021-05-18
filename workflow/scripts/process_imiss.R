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
f_list <- snakemake@input

#f_list <- list.files(path = path, pattern = "*.sample_missingness.imiss") %>%
#  paste0(path, "/", .)

message("Reading Missingness")

df <- tibble(INDV = character())
for (imiss in f_list) {
  chrom <- str_match(basename(imiss), "(chr\\d+).")[2]
  df <- read_tsv(imiss, col_types = cols(
    .default = "i",
    INDV = "c",
    F_MISS = "d")) %>%
    select(INDV, !!chrom := F_MISS) %>%
    left_join(df, by = "INDV")
}
message("Processing Missingness")

df %<>% select_if(function(col) max(col) >= threshold ||
                    str_detect(names(.), "INDV")) %>%
  mutate(maxs = apply(Filter(is.numeric, .), 1, max)) %>%
  filter(maxs >= threshold) %>%
  arrange(desc(maxs))

write_tsv(df, path = paste0(path, "/chrall.imiss"))
write_lines(df$INDV, paste0(path, "/chrall.irem"))
