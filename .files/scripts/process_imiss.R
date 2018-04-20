library(tidyverse)
library(magrittr)

args <- commandArgs(TRUE)

path <- path.expand(args[1])
threshold <- as.double(args[2])

f_list <- list.files(path = path, pattern = "*.sample_missingness.imiss")

df <- tibble(INDV = character())
for (imiss in f_list) {
  chrom = str_match(imiss, "(chr\\d+).")[2]
  df <- read_tsv(paste0(path,"/",imiss),
                          col_types = cols(
                                           .default = "i",
                                           INDV = "c",
                                           F_MISS = "d")) %>%
    select(INDV, !!chrom := F_MISS) %>%
    left_join(df, by = "INDV")
}

df %<>% select_if(function(col) max(col) >= threshold ||
                    str_detect(names(.), "INDV")) %>%
  mutate(maxs = apply(Filter(is.numeric, .), 1, max)) %>%
  filter(maxs >= threshold) %>%
  arrange(desc(maxs))

write_tsv(df, path = paste0(path, "/chrall.imiss"))
write_lines(df$INDV, paste0(path, "/chrall.irem"))
