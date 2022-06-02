#### clean and combine pre/post impute GWAS CYP4F2 variants with previously recorded dmet/luminex variants
#### this script feeds into "01_covariates.R" and "01_replication_covariates.R" files
#### heidiesteiner@email.arizona.edu
#### Jan 30 2022


#### change from VKOR ! 


#### load libraries 
library(tidyverse)
library(readxl)

############## 
#### AZ #####
##############

#### cyp4f2 rs2108622 GRCh37.p13 chr 19	15990431 C>T
####                  GRCh38.p13 chr 19	15879621 C>G 

#### rs2108622 AND exm1438764 in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr19:15879621:C:T  in phased_autosomes_lifted.bim 

#### vkor 1173 rs9934438 chr16:31093557:G:A for GRCh38 
####                  or chr16:31104878 for GRCh37 

#### 16:31104878-G-A in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr16:31093557:G:A in phased_autosomes_lifted.bim


#### plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted --snp chr19:15879621:C:T --recode --out cyp4f2_rs2108622_imputed_tucson
#### plink1.9 --bfile ../data/rawdata/data --snps exm1438764, rs2108622 --recode --out cyp4f2_rs2108622_raw_tucson

#### load raw GWAS data 
raw = read_table("results/datasets/cyp4f2_rs2108622_raw_tucson.ped",
                 col_names = F) %>% 
  select(rawid = X2,
         snp1 = X7, 
         snp2 = X8) %>% 
  unite("vip", starts_with("snp")) %>%
  mutate(vip = factor(vip, levels = c("0_0", "T_T", "T_C", "C_C"), labels = c("NA", "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         method = "gwas_raw") %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) 

table(raw$vip)

#### manually remove duplicate IDs
#### choosing the least missing row
raw_dups_ids = raw %>% filter(row > 1) %>% select(rawid)

#### choose everyone with mismatch calls 
raw_dups = raw %>% 
  filter(rawid %in% raw_dups_ids$rawid) %>% 
  group_by(rawid) %>% 
  filter(vip != "NA")

#### remove IDs w/ mismatch calls 
raw_clean = raw %>% 
  anti_join(raw_dups) %>% 
  filter(!(vip == "NA" & row == "2")) %>% 
  select(-row) %>% 
  filter(!grepl("NA", rawid)) %>% 
  mutate(method = factor(method))

rm(raw_dups_ids)
rm(raw_dups)
rm(raw)

#### load imputed GWAS data 
imputed = read_table("results/datasets/cyp4f2_rs2108622_imputed_tucson.ped",
                     col_names = F) %>% 
  select(id = X2,
         snp1 = X7, 
         snp2 = X8)%>% 
  unite("vip", starts_with("snp")) %>%
  mutate(vip = factor(vip, levels = c("T_T", "T_C", "C_C"), labels = c( "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         id2 = gsub("0_0_", "", id),
         rawid = gsub("_.*", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         vip,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(imputed$vip)


#### load IWPC data (contains taqman genotypes)
taqman = read_xlsx("data/UAZ.xlsx",
                   sheet = "Subject Data") %>% 
  mutate(rawid = toupper(`Subject ID`),
         method = "taqman",
         vip = factor(`CYP4F2 genotype: 1297G>A; Val433Met; chr19:15879621; rs2108622`,
                      levels = c("AA", "AG", "GG"),
                      labels = c("2", "1", "0"))) %>% 
  select(rawid, 
         vip,
         method) 

table(taqman$vip)

#### merge three datasets (taqman, raw, and imputed genotypes)
gene = raw_clean %>% 
  full_join(imputed) %>% 
  complete(rawid, method) %>% 
  full_join(taqman) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = vip) %>% 
  group_by(rawid) %>% 
  fill(gwas_raw, .direction = "downup") %>% 
  fill(gwas_imputed, .direction = "downup") %>% 
  fill(taqman, .direction = "downup") %>% 
  group_by(rawid) %>% 
  na_if(., "NA") %>%  
  select(-row) %>% 
  distinct()

table(gene$gwas_raw, gene$gwas_imputed)

write_delim(gene, "results/datasets/tucson_cyp4f2.tsv", delim = "\t")

rm(imputed)
rm(raw_clean)
rm(taqman)
rm(gene)



############## 
#### PR #####
##############

#### cyp4f2 rs2108622 GRCh37.p13 chr 19	15990431 C>T
####                  GRCh38.p13 chr 19	15879621 C>G 

#### rs2108622  in rawdata.bim
#### chr19:15879621:C:T  in phased_autosomes_lifted.bim 

#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/imputed/phased_autosomes_lifted --snps chr16:31096368:C:T, chr16:31093557:G:A --recode --out vkor_1639_1173_imputed_sj
#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/rawdata/rawdata --snps rs9923231, rs9934438 --recode --out vkor_1639_1173_raw_sj
#### scp steiner@karneslab.pharmacy.arizona.edu:/home/steiner/GWAS_MAY2021/candidate_genes/cyp2c9_2_3_*_tucson.ped ~/Documents/dissertation/GWAS/data/

#### load raw GWAS data 
raw = read_table("results/datasets/cyp4f2_rs2108622_raw_sj.ped",
                 col_names = F) %>% 
  select(id = X2,
         snp1 = X7, 
         snp2 = X8) %>% 
  unite("vip", starts_with("snp")) %>%
  mutate(vip = factor(vip, levels = c( "A_A", "A_G", "G_G"), labels = c( "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         method = "gwas_raw",
         rawid = sub(".*?_", "", id)) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  select(-id, -row)

table(raw$vip)

#### load imputed GWAS data 
imputed = read_table("results/datasets/cyp4f2_rs2108622_imputed_sj.ped",
                     col_names = F) %>% 
  select(id = X2,
         snp1 = X7, 
         snp2 = X8)%>% 
  unite("vip", starts_with("snp")) %>%
  mutate(vip = factor(vip, levels = c("T_T", "T_C", "C_C"), labels = c( "2", "1", "0")), 
         id2 = gsub("0_0_0_0_", "", id),
         rawid = sub(".*?_", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         vip,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(imputed$vip)

#### load IWPC data (contains dmet/luminex genotypes)
UPR_ids = read_excel("data/ID List.xlsx") %>%
  mutate(FID = paste0(rep(1:48,2), "_", `Sample_ID_GWAS array`)) %>%
  select(-...9, -...10, -...11, -IID, -Samples_ID_DMET_Plus_Files)

UPR <- read_excel("data/UPR.xlsx") %>%
  rename(`Patient ID Code` = 'Subject ID') %>%
  mutate(`Patient ID Code` = gsub("WPRAPR", "WPRA",`Patient ID Code`)) %>%
  inner_join(UPR_ids, by = "Patient ID Code") %>% 
  select(rawid = `Sample_ID_GWAS array`,
         vip = `CYP4F2 genotype: 1297G>A; Val433Met; chr19:15879621; rs2108622`
  ) %>% 
  mutate(id = sub(".*?_", "", rawid),
         method = "dmet_luminex",
         vip = factor(vip, levels = c("G/G", "A/G", "A/A"), labels = c("0", "1", "2"))) %>% 
  select(-id)

rm(UPR_ids)
table(UPR$vip)

#### merge three datasets (dmet/luminex, raw, and imputed genotypes)
vip = imputed %>% 
  complete(rawid, method) %>% 
  full_join(UPR) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = vip) %>% 
  group_by(rawid) %>% 
  fill(gwas_imputed, .direction = "downup") %>% 
  fill(dmet_luminex, .direction = "downup") %>% 
  group_by(rawid) %>% 
  select(-row)%>% 
  na_if(., "NA")%>% 
  distinct()  

table(vip$gwas_imputed, vip$dmet_luminex)

write_delim(vip, "results/datasets/sj_cyp4f2.tsv", delim = "\t")

rm(imputed)
rm(raw)
rm(UPR)
rm(vip)

