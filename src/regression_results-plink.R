#####  MANHATTAN PLOTS from assoc files 
##### H Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-02-02

#### example PLINK command to lead this script

# plink1.9 --vcf ../../data/imputed/topmed/phased_autosomes_lifted.$f.vcf.gz 
# --pheno ../plink_covariates.txt --linear --pheno-name dose --allow-no-sex 
# --remove ../../Tractor/runs/ibs_nat/plink/fail_IDs.txt --maf 0.01 --hwe 1e-6  
# --out vitk.$f_PCadjusted_maf1_hwe6_rmdoseancestryoutliers --ci 0.95 
# --covar ../plink_covariates.txt --covar-name PC1, PC2, PC3 --hide-covar 
# --extract range ../../candidate_genes/vitkgenes_plink.txt --const-fid


#### load packages
library(readr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)

#### load data desktop
#### files have the SNP name in GRCh38 and the location is in GRCh37 !!! 

dat_az <- fread("results/datasets/vitk.99.PCadjusted_maf1_hwe6_rmdoseancestryoutliers.assoc.linear") 
dat_pr <- fread("results/datasets/pr_vitk.99.PCadjusted_maf1_hwe6_rmdoseancestryoutliers.assoc.linear") 

don <- dat_az %>% 
  
  mutate(cohort = "Tucson") %>% 
  select(term = SNP, estimate = BETA, std.error = SE, p.value = P, conf.low = L95, conf.high = U95, cohort) %>% 
  full_join(dat_pr %>% mutate(cohort = "San Juan") %>% select(term = SNP, estimate = BETA, std.error = SE, p.value = P, conf.low = L95, conf.high = U95, cohort)) %>% 
  group_by(term) %>% 
  add_count() %>% 
  filter(n ==2) %>% # only keep SNPs genotyped in both cohorts 
  arrange(term) %>% 
  filter(all(p.value < 0.05)) %>% # only keep SNPs below 0.0125 in at least one cohort
  mutate(betamatch = if_else((estimate > 0 & lag(estimate)>0) |
                               (estimate < 0 & lag(estimate) <0),
                             "yes", "no")) %>% 
  fill(betamatch, .direction = "up") %>% 
  group_by(term) %>% 
  filter(all(betamatch == "yes") ) 

