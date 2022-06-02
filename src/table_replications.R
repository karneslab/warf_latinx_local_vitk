#### Create Table of replicates and IWPC adjustments
#### follows replicates.R
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-03-14

#### in table_clean, results are duplicated for different alleles in PC analysis 


#### load libraries
library(arsenal)
library(stringr)
library(tidyverse)

#### load data
hits = read_tsv("results/datasets/vitk_pcadj_replications_iwpc.tsv") %>% 
  distinct()
stats = read_tsv("results/datasets/vitk_pcadj_regression_replicates.tsv") 


tractor_hits = read_tsv("results/datasets/vitk_tractor_replications_iwpc.tsv")%>% 
  mutate(analysis = "Local Ancestry Adjusted Regressions")
tractor_stats = read_tsv("results/datasets/vitk_tractor_regression_replicates.tsv")



#### create PLINK data
dat = hits %>% 
  mutate(SNP = str_sub(term, end = -2)) %>% 
  mutate_all(as.character) %>% 
  dplyr::select(gene, allele,rsid, SNP,HWE = P, cohort, estimate, std.error,p.value, term) %>% 
  inner_join(stats  %>%  mutate(SNP = str_sub(term, end = -3),
                                term = sub(".(.)$","\\1",term)) %>% 
               mutate_all(as.character), 
             by = c("SNP", "cohort", "term"),
             suffix = c("?iwpc", "?PC")) %>% 
  pivot_longer(cols = (ends_with("?iwpc") | ends_with("?PC") ) ,
               names_to = c("stat", "analysis"),
               names_pattern = "(.*)\\?(.*)") %>% 
  distinct() %>% 
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  dplyr::select(gene, rsid, SNP, HWE, allele, cohort, analysis, estimate, std.error, p.value) 

#### create Tractor data 
tractor_dat = tractor_stats %>% 
  dplyr::select(SNP, cohort, starts_with("anc0"), starts_with("anc1"), starts_with("anc2")) %>% 
  pivot_longer(cols = starts_with("anc"),
               names_to = c("anc", 
                            "stat"),
               names_sep = "\\.") %>% 
  right_join(
    tractor_hits %>% 
      mutate(SNP = str_sub(term, end = -2),
             SNP = gsub( "\\.", ":", SNP))%>% 
      dplyr::select(term, SNP, cohort, beta = estimate, standard_error = std.error,
             p_value = p.value) %>% 
      pivot_longer(cols = c("beta", "standard_error", "p_value"),
                   names_to = "stat"), 
    suffix = c("?tractor", "?iwpc") ,
    by = c("SNP", "cohort", "stat")
    ) %>% 
  complete(term,SNP, cohort, anc, stat, explicit = T) %>% 
  pivot_longer(cols = starts_with("value"),
                names_prefix = "value\\?",
               names_to = c('analysis')) %>% 
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  dplyr::select(SNP, cohort, analysis, estimate = beta, std.error = standard_error, p.value=p_value, anc) 
  


#### create table data 
table = dat %>% 
  mutate(anc = "NA",
         ugh = "first") %>% 
  mutate_all(as.character) %>% 
  full_join(tractor_dat %>% mutate(ugh = "second", HWE = "oops", rsid = "also-oops") %>% mutate_all(as.character) %>% separate(SNP, into = c("chr", "pos", "a1", "a2"), remove = F) %>% 
              unite("allele", a1:a2)) %>% 
  mutate(anc = factor(anc, 
                      levels = c("anc0", "anc1", "anc2"),
                      labels = c("EUR", "NAT", "AFR")),
         SNP = fct_relevel(SNP,rev),
         SNP = gsub( "\\.", ":", SNP)) %>%
  mutate_at(.vars = c("estimate", "std.error", "p.value"), list( ~ as.numeric(.))) %>% 
  mutate_at(.vars = c("estimate", "std.error", "p.value"), list( ~ round(., digits = 3))) %>% 
  dplyr::select(gene, rsid, SNP, HWE, allele, cohort, analysis, anc, BETA = estimate, SE = std.error, P = p.value, ugh) 



#### clean up the table 
plus_minus <- "\U00B1"


table_clean = table %>% 
  ungroup() %>% 
  mutate(Allele = gsub("_", "", allele),
         SNP = str_sub(SNP, end = -5),
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
         )
  ) %>% 
  group_by(gene) %>% 
  dplyr::select(gene, SNP, rsid, HWE, Allele, Cohort, Model, Ancestry = anc, `Effect Size`, `P-value` = P, ugh)%>% 
  group_by(ugh,SNP, gene, rsid,Allele, Ancestry, Model, Cohort) %>% 
  arrange(.by_group = T) %>% 
  ungroup() %>% 
  distinct() 




#### write out tables

table_clean %>% 
  filter(ugh != "first") %>% 
  select(-ugh, -SNP) %>% 
  mutate(`P-value` = if_else(is.na(`P-value`), "-", as.character(`P-value`))) %>% 
  write_tsv("results/tables/secondary_replications.tsv")

table_clean %>% 
  filter(ugh == "first") %>% 
  select(-Ancestry, -ugh, -SNP) %>% 
  write_tsv("results/tables/primary_replications.tsv")

