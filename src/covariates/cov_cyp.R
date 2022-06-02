#### clean and combine pre/post impute GWAS CYP2C9*2/*3 variants with previously recorded taqman variants
#### this script feeds into "updatecovariates.R" and "replication_covar.R" files
#### heidiesteiner@email.arizona.edu
#### Jan 28 2022
#### updated Jan 31 2022

#### load libraries 
library(tidyverse)
library(readxl)

############## 
#### AZ #####
##############

#### cyp*2 rs1799853 chr10:94942290:C:T for GRCh38 
####              or chr10:96702047:C:T for GRCh37 

#### Seq-rs1799853-TOP-2 in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr10:94942290:C:T  in phased_autosomes_lifted.bim 

#### cyp*3 rs1057910 chr10:94981296:A:C for 38 
####              or chr10:96741053:A:G for 37 

#### exm844046 in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr10:94981296:A:C in phased_autosomes_lifted.bim


#### plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted --snps chr10:94942290:C:T, chr10:94981296:A:C --recode --out cyp2c9_2_3_imputed_tucson
#### plink1.9 --bfile ../data/rawdata/data --d : --snps Seq-rs1799853-TOP-2, exm844046 --recode --out cyp2c9_2_3_raw_tucson
#### scp steiner@karneslab.pharmacy.arizona.edu:/home/steiner/GWAS_MAY2021/candidate_genes/cyp2c9_2_3_*_tucson.ped ~/Documents/dissertation/GWAS/data/

#### load raw GWAS data 
cyp_raw = read_table("results/datasets/cyp2c9_2_3_raw_tucson.ped",
                      col_names = F) %>% 
  select(rawid = X2,
         cyp_2_snp1 = X7, 
         cyp_2_snp2 = X8,
         cyp_3_snp1 = X9,
         cyp_3_snp2 = X10)%>% 
  unite("cyp_2", starts_with("cyp_2")) %>%
  unite("cyp_3", starts_with("cyp_3")) %>% 
  mutate(cyp_2 = factor(cyp_2, levels = c("0_0", "C_C", "T_C", "T_T"), labels = c("NA", "0", "1", "2")), ## these are now coded as # of copies of variant alleles 
        cyp_3 = factor(cyp_3, levels = c("0_0", "A_A", "C_A", "C_C"), labels = c("NA", "0", "1", "2")),
        cyp_stars = if_else(cyp_2 == "0" & cyp_3 == "0", "*1/*1",
                    if_else(cyp_2 == "0" & cyp_3 %in% c("1","2"), "*1/*3",
                            if_else(cyp_2 %in% c("1","2"), "*1/*2", "-99"))),
        method = "gwas_raw") %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) 

table(cyp_raw$cyp_2, cyp_raw$cyp_3)
table(cyp_raw$cyp_stars)


#### manually remove duplicate IDs
#### choosing the least missing row
#### WARFER031 has two different cyp_2 genotypes 
raw_dups_ids = cyp_raw %>% filter(row > 1) %>% select(rawid)
warfer031 = cyp_raw %>% filter(rawid == "WARFER031" & cyp_stars == "*1/*1")

cyp_raw_dups = cyp_raw %>% 
  filter(rawid %in% raw_dups_ids$rawid) %>% 
  filter(cyp_stars == "-99") %>% 
  full_join(warfer031)

cyp_raw_clean = cyp_raw %>% 
  anti_join(cyp_raw_dups) %>% 
  select(-row) %>% 
  filter(!grepl("NA", rawid)) %>% 
  mutate(method = factor(method))

rm(raw_dups_ids)
rm(warfer031)
rm(cyp_raw_dups)
rm(cyp_raw)


#### load imputed GWAS data 
cyp_imputed = read_table("results/datasets/cyp2c9_2_3_imputed_tucson.ped",
                     col_names = F) %>% 
  select(id = X2,
         cyp_2_snp1 = X7, 
         cyp_2_snp2 = X8,
         cyp_3_snp1 = X9,
         cyp_3_snp2 = X10)%>% 
  unite("cyp_2", starts_with("cyp_2")) %>%
  unite("cyp_3", starts_with("cyp_3")) %>% 
  mutate(cyp_2 = factor(cyp_2, levels = c("C_C", "T_C", "T_T"), labels = c( "0", "1", "2")), ## these are now coded as # of copies of variant alleles 
         cyp_3 = factor(cyp_3, levels = c( "A_A", "C_A", "C_C"), labels = c( "0", "1", "2")),
         cyp_stars = if_else(cyp_2 == "0" & cyp_3 == "0", "*1/*1",
                             if_else(cyp_2 == "0" & cyp_3 %in% c("1","2"), "*1/*3",
                                     if_else(cyp_2 %in% c("1","2"), "*1/*2", "-99"))),
         id2 = gsub("0_0_", "", id),
         rawid = gsub("_.*", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         cyp_2,
         cyp_3,
         cyp_stars,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(cyp_imputed$cyp_2, cyp_imputed$cyp_3)
table(cyp_imputed$cyp_stars)

#### load IWPC data (contains taqman genotypes)
cyp_taqman = read_xlsx("data/UAZ.xlsx",
                      sheet = "Subject Data") %>% 
  mutate(rawid = toupper(`Subject ID`),
         method = "taqman",
         cyp_2 = NA,
         cyp_3 = NA) %>% 
  select(rawid, 
         cyp_2,
         cyp_3,
         cyp_stars = `Cyp2C9 genotypes`,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method =factor(method)) %>% 
  select(-row)

table(cyp_taqman$cyp_stars)

#### merge three datasets (taqman, raw, and imputed genotypes)
cyps = cyp_raw_clean %>% 
  full_join(cyp_imputed) %>% 
  complete(rawid, method) %>% 
  full_join(cyp_taqman) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
             values_from = c("cyp_2", "cyp_3", "cyp_stars")) %>% 
  group_by(rawid) %>% 
  fill(cyp_2_gwas_raw, .direction = "downup") %>% 
  fill(cyp_2_gwas_imputed, .direction = "downup") %>% 
  fill(cyp_2_taqman, .direction = "downup") %>% 
  fill(cyp_3_gwas_raw, .direction = "downup") %>% 
  fill(cyp_3_gwas_imputed, .direction = "downup") %>% 
  fill(cyp_3_taqman, .direction = "downup") %>% 
  fill(cyp_stars_gwas_raw, .direction = "downup") %>% 
  fill(cyp_stars_gwas_imputed, .direction = "downup")%>% 
  fill(cyp_stars_taqman, .direction = "downup") %>% 
  group_by(rawid)  %>% 
  na_if(., "NA") %>% 
  mutate(cyp_stars_gwas_raw = na_if(cyp_stars_gwas_raw, "-99"))%>% 
  select(-row) %>% 
  distinct() %>% 
  discard(~all(is.na(.))) 

table(cyps$cyp_stars_gwas_imputed, cyps$cyp_stars_gwas_raw, cyps$cyp_stars_taqman)
table(cyps$cyp_2_gwas_raw, cyps$cyp_2_gwas_imputed)
table(cyps$cyp_3_gwas_raw, cyps$cyp_3_gwas_imputed)

write_delim(cyps, "results/datasets/tucson_cyp.tsv", delim = "\t")

rm(cyp_imputed)
rm(cyp_raw_clean)
rm(cyp_taqman)
rm(cyps)


############## 
#### PR #####
##############

#### cyp*2 rs1799853 chr10:94942290:C:T for GRCh38 
####              or chr10:96702047:C:T for GRCh37 

#### not present in rawdata.bim
#### chr10:94942290:C:T  in phased_autosomes_lifted.bim 

#### cyp*3 rs1057910 chr10:94981296:A:C for 38 
####              or chr10:96741053:A:G for 37 

#### not present in rawdata.bim
#### chr10:94981296:A:C in phased_autosomes_lifted.bim


#### plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted --snps chr10:94942290:C:T, chr10:94981296:A:C --recode --out cyp2c9_2_3_imputed_tucson
#### scp steiner@karneslab.pharmacy.arizona.edu:/home/steiner/GWAS_MAY2021/candidate_genes/cyp2c9_2_3_*_sj.ped ~/Documents/dissertation/GWAS/data/

#### load imputed GWAS data 
cyp_imputed = read_table("results/datasets/cyp2c9_2_3_imputed_sj.ped",
                         col_names = F) %>% 
  select(id = X2,
         cyp_2_snp1 = X7, 
         cyp_2_snp2 = X8,
         cyp_3_snp1 = X9,
         cyp_3_snp2 = X10)%>% 
  unite("cyp_2", starts_with("cyp_2")) %>%
  unite("cyp_3", starts_with("cyp_3")) %>% 
  mutate(cyp_2 = factor(cyp_2, levels = c("C_C", "T_C", "T_T"), labels = c( "0", "1", "2")), ## these are now coded as # of copies of variant alleles 
         cyp_3 = factor(cyp_3, levels = c( "A_A", "C_A", "C_C"), labels = c( "0", "1", "2")),
         cyp_stars = if_else(cyp_2 == "0" & cyp_3 == "0", "*1/*1",
                             if_else(cyp_2 == "0" & cyp_3 %in% c("1","2"), "*1/*3",
                                     if_else(cyp_2 %in% c("1","2"), "*1/*2", "-99"))),
         id2 = gsub("0_0_0_0_*_", "", id),
         rawid = sub(".*?_", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         cyp_2,
         cyp_3,
         cyp_stars,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(cyp_imputed$cyp_2, cyp_imputed$cyp_3)
table(cyp_imputed$cyp_stars)

#### load IWPC data (contains dmet/luminex genotypes)
UPR_ids = read_excel("data/ID List.xlsx") %>%
  mutate(FID = paste0(rep(1:48,2), "_", `Sample_ID_GWAS array`)) %>%
  select(-...9, -...10, -...11, -IID, -Samples_ID_DMET_Plus_Files)

UPR <- read_excel("data/UPR.xlsx") %>%
  rename(`Patient ID Code` = 'Subject ID') %>%
  mutate(`Patient ID Code` = gsub("WPRAPR", "WPRA",`Patient ID Code`)) %>%
  inner_join(UPR_ids, by = "Patient ID Code") %>% 
  select(rawid = `Sample_ID_GWAS array`,
         cyp_stars = `Cyp2C9 genotypes`) %>% 
  mutate(
         method = "dmet_luminex")

rm(UPR_ids)

#### merge two datasets (imputed and dmet/luminex)
cyps_pr = cyp_imputed %>% 
  full_join(UPR) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = c("cyp_2", "cyp_3", "cyp_stars")) %>% 
  group_by(rawid) %>% 
  fill(cyp_2_gwas_imputed, .direction = "downup") %>% 
  fill(cyp_3_gwas_imputed, .direction = "downup") %>% 
  fill(cyp_stars_gwas_imputed, .direction = "downup")%>% 
  fill(cyp_stars_dmet_luminex, .direction = "downup") %>% 
  group_by(rawid)  %>% 
  na_if(., "NA") %>% 
  mutate(cyp_stars_dmet_luminex_23only = if_else(cyp_stars_dmet_luminex == "*1/*1", "*1/*1", 
                                 if_else(cyp_stars_dmet_luminex %in% c("*1/*11", "*1/*2", "*1/*3",
                                                                       "*1/*4", "*1/*8", "*1/*9"), "*1/*decrease",
                                         if_else(cyp_stars_dmet_luminex %in% c("*2/*3", "*2/*9", 
                                                                               "*3/*5", "*3/*6", "*3/*8"), "*decrease/*decrease", "-99"))))%>% 
  select(-row) %>% 
  distinct() %>% 
  discard(~all(is.na(.))) 

rm(cyp_taqman)
rm(cyp_imputed)
rm(UPR)

table(cyps_pr$cyp_stars_gwas_imputed,  cyps_pr$cyp_stars_dmet_luminex)
table(cyps_pr$cyp_stars_gwas_imputed,  cyps_pr$cyp_stars_dmet_luminex_23only)

write_delim(cyps_pr, "results/datasets/sj_cyp.tsv", delim = "\t")


rm(cyps_pr)
