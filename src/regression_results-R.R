#### Create Table of results from PC adjusted linear models
#### follows lm.R
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Created 2022-03-12

#### load libraries
library(tidyverse)

#### load data
lm2 = read_tsv("results/datasets/vitk_unadj.tsv") %>% 
  mutate(SNP = str_sub(term, end = -3),
         term = as.factor(term)) %>% 
  complete(term) %>% 
  group_by(term) %>% 
  add_count() %>% 
  filter(n ==2) %>% # only keep SNPs genotyped in both cohorts 
  arrange(term) %>% 
  filter(any(p.value < 0.0125 & lag(p.value)< 0.0125))%>% # only keep SNPs below 0.0125 in both cohorts
  mutate(betamatch = if_else((estimate > 0 & lag(estimate)>0) |
                               (estimate < 0 & lag(estimate) <0),
                             "yes", "no"))%>% 
  fill(betamatch, .direction = "up") %>% 
  group_by(SNP) %>% 
  filter(all(betamatch == "yes") ) 


length(unique(lm2$SNP)) # SNPs that "replicate" at 0.0125 

######## plink query 
lm2  %>% 
  ungroup() %>% 
  dplyr::select(SNP) %>% 
  mutate(SNP = gsub("\\.", "\\:", SNP)) %>% 
  unique() %>% 
  write_tsv(., "results/datasets/vitk_pcadj_replications.txt")

#### write out summary stats of replicated variants for table_replications.R
#### this dataset features all the variants with p < 0.0125 in BOTH cohorts and matching direction of effect (beta) 
lm2 %>% 
  select(term, estimate, std.error, p.value, cohort) %>% 
  write_tsv("results/datasets/vitk_pcadj_regression_replicates.tsv")  



