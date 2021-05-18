library(tidyverse)

demo.path = snakemake@input[['demo']]
fam.path = snakemake@input[['fam']]
out.path = snakemake@output[['out']]

demo.raw <- read_csv(demo.path)
fam.raw <- read_table2(fam.path, col_names = F)

demo <- demo.raw %>%
  select(RID, PTID, VISCODE, EXAMDATE, DX, AGE, PTGENDER, MMSE) %>%
  arrange(RID, EXAMDATE) %>%
  group_by(PTID) %>%
  filter(!is.na(MMSE)) %>%
  slice(which.max(EXAMDATE)) %>%
  ungroup() %>%
  mutate(sex = recode(PTGENDER, Male = 1, Female = 2)) %>%
  select(FID = PTID, IID = PTID, sex, MMSE)


fam <- fam.raw %>%
  rename("FID" = X1, "IID" = X2, "MIID" = X3, "PIID" = X4, "sex" = X5, "pheno" = X6) %>%
  mutate(FID = str_extract(FID, "(?<=_)[:alnum:].*"),
         IID = str_extract(IID, "(?<=_)[:alnum:].*")) %>%
  select(-sex, -pheno) %>%
  left_join(., demo, by = c("FID", "IID")) %>%
  select(FID, IID, PIID, MIID, sex, MMSE) %>%
  mutate(sex = replace_na(sex, 0),
         sex = as_factor(sex),
         MMSE = replace_na(MMSE, -9)) %>%
 mutate(MMSE = -9)

write_delim(fam, out.path, col_names = F)
