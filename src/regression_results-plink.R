#### IN: linear.assoc files direct from PLINK
#### OUT: plink-query.txt and .tsv input to table_replications.R
#### H Steiner
#### heidiesteiner@email.arizona.edu
#### Created: 2022-02-02
#### Updated: 2022-06-19

#### example PLINK command to input into this script

# for f in 80 90 99; do 
# plink1.9 --vcf ../../data/imputed/topmed/phased_autosomes_lifted.$f.vcf.gz 
# --pheno ../plink_covariates.txt --linear hide-covar --pheno-name dose 
# --allow-no-sex --remove ../../Tractor/runs/ibs_nat/plink/fail_IDs.txt 
# --maf 0.05 --hwe 1e-6 --out vitk.$f.PCadjusted_maf5_hwe6_rmdoseoutliers 
# --ci 0.95 --covar ../plink_covariates.txt --covar-name PC1, PC2, PC3 
# --extract range ../../candidate_genes/vitkgenes_plink.txt --const-fid
# ;done


#### load packages
library(readr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)

#### load data desktop
#### files have the SNP name in GRCh38 and the location is in GRCh37 !!! 

dat_az <- fread("results/datasets/vitk.90.PCadjusted_maf5_hwe6_rmdoseoutliers.assoc.linear") 
dat_pr <- fread("results/datasets/pr_vitk.90.PCadjusted_maf5_hwe6.assoc.linear") 

don <- dat_az %>% 
  
  mutate(cohort = "Tucson") %>% 
  select(term = SNP, estimate = BETA, std.error = SE, p.value = P, conf.low = L95, conf.high = U95, cohort) %>% 
  full_join(dat_pr %>% mutate(cohort = "San Juan") %>% select(term = SNP, estimate = BETA, std.error = SE, p.value = P, conf.low = L95, conf.high = U95, cohort)) %>% 
  group_by(term) %>% 
  add_count() %>% 
  filter(n ==2) %>% # only keep SNPs genotyped in both cohorts 
  arrange(term) %>% 
  filter(all(p.value < 0.0125)) %>% # only keep SNPs below 0.0125 in both cohorts
  mutate(betamatch = if_else((estimate > 0 & lag(estimate)>0) |
                               (estimate < 0 & lag(estimate) <0),
                             "yes", "no")) %>% 
  fill(betamatch, .direction = "up") %>% 
  group_by(term) %>% 
  filter(all(betamatch == "yes") ) 


#### write out a plink query for testing replications
don  %>% 
  ungroup() %>% 
  dplyr::select(SNP = term) %>% 
  unique() %>% 
  write_tsv(., "results/datasets/vitk90_pcadj_maf5_replications.txt")

#### write out summary stats of replicated variants for table_replications.R
#### this dataset features all the variants with p < 0.0125 in BOTH cohorts and matching direction of effect (beta) 
don %>% 
  select(term, estimate, std.error, p.value, cohort) %>% 
  write_tsv("results/datasets/vitk90_pcadj_maf5_regression_replicates.tsv")  




