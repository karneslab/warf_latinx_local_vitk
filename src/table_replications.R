#### IN: replications_iwpc.tsvs from replicates.R and tractor_replicates.R
#### IN: replications.tsvs from regression-results.R and tractor_regression-results.R
#### OUT: Tables for pasting in manuscript (almost complete)
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Created: 2022-03-14
#### Updated: 2022-06-22




#### load libraries
library(arsenal)
library(stringr)
library(tidyverse)

options(pillar.sigfig = 3)

#### load plink data
hits = read_tsv("results/datasets/vitk90_maf05_replications_iwpcadj.tsv") %>% ## after IWPC adjustment
  distinct()
stats = read_tsv("results/datasets/vitk90_pcadj_maf5_regression_replicates.tsv") ## results from regression_results.R


### load tractor data 
tractor_hits = read_tsv("results/datasets/vitk_tractor_replications_iwpc.tsv")%>% ## after IWPC adjustment
  mutate(analysis = "Local Ancestry Adjusted Regressions")
tractor_stats = read_tsv("results/datasets/vitk90_LAadj_maf5_regression_replicates.tsv")



#### create PLINK data
dat = hits %>% 
  mutate(SNP = str_sub(term, end = -1),
         allele = gsub("_", "", allele)) %>% 
  unite("chr_pos", chr:pos) %>% 
  # mutate_all(as.character) %>% 
  dplyr::select(chr_pos, allele, rsid, SNP,HWE = P, cohort, estimate, std.error,p.value, term) %>% 
  inner_join(stats  %>%  mutate(SNP = str_sub(term, end = -1),
                                term = sub(".(.)$","\\1",term)) %>% 
               separate(term, into = c("chr", "pos", "allele"), sep = ":") %>% 
               unite( "chr_pos", chr:pos, remove = T, sep = "_") , 
             by = c("allele", "cohort", "chr_pos"),
             suffix = c("?iwpc", "?PC")) %>% 
  dplyr::select(-starts_with("SNP")) %>% 
  pivot_longer(cols = (ends_with("?iwpc") | ends_with("?PC") ) ,
               names_to = c("stat", "analysis"),
               names_pattern = "(.*)\\?(.*)") %>% 
  distinct() %>% 
  pivot_wider(names_from = stat,
              values_from = value) 

#### create Tractor data 
tractor_dat = tractor_stats %>% 
  dplyr::select(SNP, cohort, starts_with("anc0"), starts_with("anc1"), starts_with("anc2")) %>% 
  pivot_longer(cols = starts_with("anc"),
               names_to = c("anc", 
                            "stat"),
               names_sep = "\\.") %>% 
 # mutate(stat = factor(stat, levels = c("beta", "p_value", "standard_error", "t_stat", "HWE"))) %>% 
  right_join(
    tractor_hits %>% 
      mutate(SNP = str_sub(term, end = -1),
             SNP = gsub( "\\.", ":", SNP),
             anc = "-")%>% 
      dplyr::select( SNP, cohort, beta = estimate, standard_error = std.error,
             p_value = p.value, HWE = P) %>% 
      pivot_longer(cols = c("beta", "standard_error", "p_value"),
                   names_to = "stat"), 
    #  mutate(stat = factor(stat, levels = c("beta", "p_value", "standard_error", "t_stat", "HWE"))), 
    suffix = c("?tractor", "?iwpc") ,
    by = c("SNP", "cohort", "stat")
    ) %>% 
  pivot_longer(cols = starts_with("value"),
                names_prefix = "value\\?",
               names_to = c('analysis')) %>% 
  distinct() %>% 
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  distinct(across(c(SNP,beta)), .keep_all = T) %>% 
  dplyr::select(chr_pos = SNP, cohort, analysis, estimate = beta, std.error = standard_error, p.value=p_value, anc, HWE) 
  
view = tractor_dat %>% 
  filter(analysis == "iwpc", p.value < 0.0125)

view %>% 
  group_by(cohort) %>% 
  summarise(max = max(abs(estimate)),
            min = min(abs(estimate)))

View(view)
#### create table data 
table = dat %>% 
  mutate(anc = "NA",
         ugh = "first") %>% 
  mutate_all(as.character) %>% 
  full_join(tractor_dat %>% 
              mutate(ugh = "second", rsid = "also-oops") %>% 
              mutate_all(as.character) %>% 
              separate(chr_pos, into = c("chr", "pos", "a1", "a2"), remove = T) %>% 
              unite("allele", a1:a2) %>% 
              unite("chr_pos", chr:pos, sep = "_")) %>% 
  mutate(anc = factor(anc, 
                      levels = c("anc0", "anc1", "anc2"),
                      labels = c("EUR", "NAT", "AFR")),
         chr_pos = fct_relevel(chr_pos,rev),
         chr_pos = gsub( "\\.", ":", chr_pos)) %>%
  mutate_at(.vars = c("estimate", "std.error", "p.value"), list( ~ as.numeric(.))) %>% 
  mutate_at(.vars = c("estimate", "std.error"), list( ~ round(., digits = 2))) %>% 
  mutate(p.value = formatC(p.value, format = "e", digits = 2),
         p.value = gsub("e","x10", p.value),
         p.value = gsub("-0", "-", p.value)) %>% 
  dplyr::select(rsid, SNP = chr_pos, HWE, allele, cohort, analysis, anc, BETA = estimate, SE = std.error, P = p.value, ugh) 

View(table)

#### clean up the table 
plus_minus <- "\U00B1"



table_clean = table %>% 
  ungroup() %>% 
  mutate(Allele = gsub("_", "", allele),
         # SNP = str_sub(SNP, end = -4),
         SNP = gsub(":$", "", SNP),
         SNP = fct_relevel(SNP, rev),
         Cohort = fct_relevel(cohort, rev),
         Model = factor(analysis, 
                        levels = c("PC",
                                   "tractor",
                                   "iwpc"),
                        labels = c("PC",
                                   "Tractor",
                                   "IWPC")),
         `Effect Size` = paste(BETA,paste0("(",plus_minus,SE,")"),
                               sep = " "),
         `Effect Size` = case_when(
           `Effect Size` == "NA (±NA)" ~ "-",
           TRUE ~ `Effect Size`
         ),
         anc = if_else(Model == "IWPC" & ugh == "second", "-", as.character(anc)),
         gene = if_else(str_starts(SNP, "chr1_656"), "AK4", 
                        if_else(str_starts(SNP,  "chr2_750"), "HK2",
                                if_else(str_starts(SNP, "chr5_180"), "MGAT1",
                                        if_else(str_starts(SNP, "chr16_31"), "VKORC1", 
                                            if_else(str_starts(SNP, "chr10_95"), "ACSM6, LOC107984257", 
                                                    if_else(str_starts(SNP, "chr17_1420"), "COX10", 
                                                            if_else(str_starts(SNP, "chr11.747"), "NEU3/OR2AT2P",
                                                                    if_else(str_starts(SNP, "chr18.442"), "ST8SIA5", "No gene identified"))))))))
  ) %>% 
  group_by(gene) %>% 
  dplyr::select(gene, rsid, SNP, HWE, Allele, Cohort, Model, Ancestry = anc, `Effect Size`, `P-value` = P, ugh)%>% 
  group_by(gene, ugh,SNP,Allele, Ancestry, Model, Cohort) %>% 
  arrange(.by_group = T) %>% 
  ungroup() %>% 
  distinct() 



View(table_clean)

#### write out tables

table_clean %>% 
  filter(ugh != "first") %>% 
  select(-ugh, -rsid, -Allele) %>% 
  mutate(`P-value` = if_else(is.na(`P-value`), "-", as.character(`P-value`))) %>% 
  group_by(gene,SNP,Model, Ancestry, Cohort) %>% 
  arrange(.by_group = T) %>% 
  write_tsv("results/tables/secondary_replications.tsv") 



table_clean %>% 
  filter(ugh == "first") %>% 
  unite("new_SNP", SNP, rsid, sep = " (") %>% 
  mutate(SNP = paste0(new_SNP, ")")) %>% 
  select(gene, SNP, Allele, HWE, Cohort, Model, `Effect Size`, `P-value`) %>% 
  write_tsv("results/tables/primary_replications.tsv")

