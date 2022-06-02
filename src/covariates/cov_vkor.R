#### clean and combine pre/post impute GWAS VKORC1 variants with previously recorded taqman/dmet/luminex variants
#### this script feeds into ")1_covariates.R" and "01_replication_covariates.R" files
#### heidiesteiner@email.arizona.edu
#### Jan 28 2022

#### consider vkor 3730 rs7294 


#### load libraries 
library(tidyverse)
library(readxl)

############## 
#### AZ #####
##############

#### vkor 1639 rs9923231 chr16:31096368:C:T for GRCh38 
####                  or chr16:31107689 for GRCh37 

#### not present in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr16:31096368:C:T  in phased_autosomes_lifted.bim 

#### vkor 1173 rs9934438 chr16:31093557:G:A for GRCh38 
####                  or chr16:31104878 for GRCh37 

#### 16:31104878-G-A in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr16:31093557:G:A in phased_autosomes_lifted.bim


#### plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted --snps chr16:31096368:C:T, chr16:31093557:G:A --recode --out vkor_1639_1173_imputed_tucson
#### plink1.9 --bfile ../data/rawdata/data --d : --snp 16:31104878-G-A --recode --out vkor_1173_raw_tucson

#### load raw GWAS data 
vkor_raw = read_table("results/datasets/vkor_1173_raw_tucson.ped",
                      col_names = F) %>% 
  select(rawid = X2,
         vkor_1173_snp1 = X7, 
         vkor_1173_snp2 = X8) %>% 
  unite("vkor_1173", starts_with("vkor_1173")) %>%
  mutate(vkor_1173 = factor(vkor_1173, levels = c("0_0", "A_A", "A_G", "G_G"), labels = c("NA", "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         method = "gwas_raw") %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) 

table(vkor_raw$vkor_1173)

#### manually remove duplicate IDs
#### choosing the least missing row

raw_dups_ids = vkor_raw %>% filter(row > 1) %>% select(rawid)
warfer007 = vkor_raw %>% filter(rawid == "WARFER007" & vkor_1173 == "NA")

vkor_raw_dups = vkor_raw %>% 
  filter(rawid %in% raw_dups_ids$rawid,
         rawid != "WARFER007") %>% 
  full_join(warfer007)

vkor_raw_clean = vkor_raw %>% 
  anti_join(vkor_raw_dups) %>% 
  select(-row) %>% 
  filter(!grepl("NA", rawid)) %>% 
  mutate(method = factor(method),
         vkor_1639 = NA)

rm(raw_dups_ids)
rm(vkor_raw_dups)
rm(vkor_raw)
rm(warfer007)

#### load imputed GWAS data 
vkor_imputed = read_table("results/datasets/vkor_1639_1173_imputed_tucson.ped",
                          col_names = F) %>% 
  select(id = X2,
         vkor_1639_snp1 = X7, 
         vkor_1639_snp2 = X8,
         vkor_1173_snp1 = X9,
         vkor_1173_snp2 = X10)%>% 
  unite("vkor_1639", starts_with("vkor_1639")) %>%
  unite("vkor_1173", starts_with("vkor_1173")) %>% 
  mutate(vkor_1639 = factor(vkor_1639, levels = c("A_A", "A_G", "G_G"), labels = c( "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         vkor_1173 = factor(vkor_1173, levels = c( "C_C", "T_C", "T_T"), labels = c( "0", "1", "2")),
         id2 = gsub("0_0_", "", id),
         rawid = gsub("_.*", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         vkor_1639,
         vkor_1173,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(vkor_imputed$vkor_1173, vkor_imputed$vkor_1639)

#### load IWPC data (contains taqman genotypes)
vkor_taqman = read_xlsx("data/UAZ.xlsx",
                       sheet = "Subject Data") %>% 
  mutate(rawid = toupper(`Subject ID`),
         method = "taqman") %>% 
  select(rawid, 
         vkor_1173 = `VKORC1 genotype:   1173 C>T (6484); chr16:31012379(hg18); rs9934438`,
         vkor_1639 = `VKORC1 genotype:   -1639 G>A (3673); chr16:31015190(hg18); rs9923231`,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method),
         vkor_1639 = factor(vkor_1639, levels = c("AA", "AG", "GG"), labels = c("2", "1", "0"))) %>% 
  select(-row)

table(vkor_taqman$vkor_1639)

#### merge three datasets (taqman, raw, and imputed genotypes)
vkors = vkor_raw_clean %>% 
  full_join(vkor_imputed) %>% 
  complete(rawid, method) %>% 
  full_join(vkor_taqman) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = c("vkor_1639", "vkor_1173")) %>% 
  group_by(rawid) %>% 
  fill(vkor_1639_gwas_raw, .direction = "downup") %>% 
  fill(vkor_1639_gwas_imputed, .direction = "downup") %>% 
  fill(vkor_1639_taqman, .direction = "downup") %>% 
  fill(vkor_1173_gwas_raw, .direction = "downup") %>% 
  fill(vkor_1173_gwas_imputed, .direction = "downup") %>% 
  fill(vkor_1173_taqman, .direction = "downup") %>% 
  select(-row)%>% 
  distinct() %>% 
  na_if(., "NA") %>% 
  discard(~all(is.na(.))) 

table(vkors$vkor_1173_gwas_raw, vkors$vkor_1173_gwas_imputed)
table( vkors$vkor_1639_gwas_imputed, vkors$vkor_1639_taqman)

write_delim(vkors, "results/datasets/tucson_vkor.tsv", delim = "\t")

rm(vkor_imputed)
rm(vkor_raw_clean)
rm(vkor_taqman)
rm(vkors)


############## 
#### PR #####
##############

#### vkor 1639 rs9923231 chr16:31096368:C:T for GRCh38 
####                  or chr16:31107689 for GRCh37 

#### rs9923231 in rawdata.bim
#### chr16:31096368:C:T  in phased_autosomes_lifted.bim 

#### vkor 1173 rs9934438 chr16:31093557:G:A for GRCh38 
####                  or chr16:31104878 for GRCh37 

#### rs9934438 in rawdata.bim
#### chr16:31093557:G:A in phased_autosomes_lifted.bim

#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/imputed/phased_autosomes_lifted --snps chr16:31096368:C:T, chr16:31093557:G:A --recode --out vkor_1639_1173_imputed_sj
#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/rawdata/rawdata --snps rs9923231, rs9934438 --recode --out vkor_1639_1173_raw_sj
#### scp steiner@karneslab.pharmacy.arizona.edu:/home/steiner/GWAS_MAY2021/candidate_genes/cyp2c9_2_3_*_tucson.ped ~/Documents/dissertation/GWAS/data/

#### load raw GWAS data 
vkor_raw = read_table("results/datasets/vkor_1639_1173_raw_sj.ped",
                      col_names = F) %>% 
  select(id = X2,
         vkor_1639_snp1 = X7, 
         vkor_1639_snp2 = X8,
         vkor_1173_snp1 = X9,
         vkor_1173_snp2 = X10) %>% 
  unite("vkor_1173", starts_with("vkor_1173")) %>%
  unite("vkor_1639", starts_with("vkor_1639")) %>% 
  mutate(vkor_1173 = factor(vkor_1173, levels = c("0_0", "A_A", "A_G", "G_G"), labels = c("NA", "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         vkor_1639 = factor(vkor_1639, levels = c("A_A", "A_G", "G_G"), labels = c( "2", "1", "0")),
         method = "gwas_raw",
         rawid = sub(".*?_", "", id)) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  select(-id, -row)

table(vkor_raw$vkor_1173)
table(vkor_raw$vkor_1639)


#### load imputed GWAS data 
vkor_imputed = read_table("results/datasets/vkor_1639_1173_imputed_sj.ped",
                          col_names = F) %>% 
  select(id = X2,
         vkor_1639_snp1 = X7, 
         vkor_1639_snp2 = X8,
         vkor_1173_snp1 = X9,
         vkor_1173_snp2 = X10)%>% 
  unite("vkor_1639", starts_with("vkor_1639")) %>%
  unite("vkor_1173", starts_with("vkor_1173")) %>% 
  mutate(vkor_1639 = factor(vkor_1639, levels = c("A_A", "A_G", "G_G"), labels = c( "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         vkor_1173 = factor(vkor_1173, levels = c( "C_C", "T_C", "T_T"), labels = c( "0", "1", "2")),
         id2 = gsub("0_0_0_0_", "", id),
         rawid = sub(".*?_", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         vkor_1639,
         vkor_1173,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(vkor_imputed$vkor_1173, vkor_imputed$vkor_1639)

#### load IWPC data (contains dmet/luminex genotypes)
UPR_ids = read_excel("data/ID List.xlsx") %>%
  mutate(FID = paste0(rep(1:48,2), "_", `Sample_ID_GWAS array`)) %>%
  select(-...9, -...10, -...11, -IID, -Samples_ID_DMET_Plus_Files)

UPR <- read_excel("data/UPR.xlsx") %>%
  rename(`Patient ID Code` = 'Subject ID') %>%
  mutate(`Patient ID Code` = gsub("WPRAPR", "WPRA",`Patient ID Code`)) %>%
  inner_join(UPR_ids, by = "Patient ID Code") %>% 
  select(rawid = `Sample_ID_GWAS array`,
         vkor_1173 = `VKORC1 genotype:   1173 C>T (6484); chr16:31012379(hg18); rs9934438`,
         vkor_1639 = `VKORC1 genotype:   -1639 G>A (3673); chr16:31015190(hg18); rs9923231`) %>% 
  mutate(
         method = "dmet_luminex",
         vkor_1173 = factor(vkor_1173, levels = c("C/C", "C/T", "T/T", "NA"), labels = c("0", "1", "2", "NA")),
         vkor_1639 = factor(vkor_1639, levels = c("A/A", "G/A",  "G/G", "NA"), labels = c("2", "1", "0", "NA")))

rm(UPR_ids)
table(UPR$vkor_1173, UPR$vkor_1639)





#### merge three datasets (dmet/luminex, raw, and imputed genotypes)
vkors = vkor_raw %>% 
  full_join(vkor_imputed) %>% 
  complete(rawid, method) %>% 
  full_join(UPR) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = c("vkor_1639", "vkor_1173")) %>% 
  group_by(rawid) %>% 
  fill(vkor_1639_gwas_raw, .direction = "updown") %>% 
  fill(vkor_1639_gwas_imputed, .direction = "downup") %>% 
  fill(vkor_1639_dmet_luminex, .direction = "downup") %>% 
  fill(vkor_1173_gwas_raw, .direction = "downup") %>% 
  fill(vkor_1173_gwas_imputed, .direction = "downup") %>% 
  fill(vkor_1173_dmet_luminex, .direction = "downup") %>% 
  group_by(rawid) %>% 
  select(-row) %>% 
  distinct() %>% 
  na_if(., "NA")

table(vkors$vkor_1173_gwas_raw, vkors$vkor_1173_gwas_imputed, vkors$vkor_1173_dmet_luminex)
table(vkors$vkor_1639_gwas_imputed, vkors$vkor_1639_gwas_raw, vkors$vkor_1639_dmet_luminex)

 write_delim(vkors, "results/datasets/sj_vkor.tsv", delim = "\t")

rm(vkor_imputed)
rm(vkor_raw)
rm(UPR)

rm(vkors)
