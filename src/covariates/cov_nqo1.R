#### clean and combine pre/post impute GWAS NQO1 variants with previously recorded taqman/dmet/luminex variants
#### this script feeds into "replication_covar.R" and "replication_covar.R" files
#### heidiesteiner@email.arizona.edu
#### Jan 30 2022



#### load libraries 
library(tidyverse)
library(readxl)

############## 
#### AZ #####
##############

#### NQO1*2 (rs) 
#### GRCh37.p13 chr 16	69745145 G>A
#### GRCh38.p13 chr 16	69711242 G>A

#### exm1253081 in deliverable_JK-5028-Karnes_BestPerformingSNPs.bim
#### chr16:69711242:G:A  in phased_autosomes_lifted.bim 

#### load raw GWAS data 
raw = read_table("results/datasets/nqo1_rs1800566_raw_tucson.ped",
                      col_names = F) %>% 
  select(rawid = X2,
         snp1 = X7, 
         snp2 = X8) %>% 
  unite("nqo1_rs1800566", starts_with("snp")) %>%
  mutate(nqo1_rs1800566 = factor(nqo1_rs1800566, levels = c("0_0", "A_A", "A_G", "G_G"), labels = c("NA", "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         method = "gwas_raw") %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) 

table(raw$nqo1_rs1800566)

#### manually remove duplicate IDs
#### choosing the least missing row
raw_dups_ids = raw %>% filter(row > 1) %>% select(rawid)

#### choose everyone with mismatch calls 
raw_dups = raw %>% 
  filter(rawid %in% raw_dups_ids$rawid) %>% 
  group_by(rawid) %>% 
  filter(nqo1_rs1800566 != "NA")

raw_clean = raw %>% 
  anti_join(raw_dups) %>% 
  filter(!(nqo1_rs1800566 == "NA" & row == "2")) %>% 
  select(-row) %>% 
  filter(!grepl("NA", rawid)) %>% 
  mutate(method = factor(method))

rm(raw_dups_ids)
rm(raw_dups)
rm(raw)

#### load imputed GWAS data 
imputed = read_table("results/datasets/nqo1_rs1800566_imputed_tucson.ped",
                          col_names = F) %>% 
  select(id = X2,
         snp1 = X7, 
         snp2 = X8)%>% 
  unite("nqo1_rs1800566", starts_with("snp")) %>%
  mutate(nqo1_rs1800566 = factor(nqo1_rs1800566, levels = c("A_A", "A_G", "G_G"), labels = c( "2", "1", "0")), ## these are now coded as # of copies of variant alleles 
         id2 = gsub("0_0_", "", id),
         rawid = gsub("_.*", "", id2),
         method = "gwas_imputed") %>% 
  select(rawid, 
         nqo1_rs1800566,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method)) %>% 
  select(-row)

table(imputed$nqo1_rs1800566)

#### load IWPC data (contains taqman genotypes)
taqman = read_xlsx("data/UAZ.xlsx",
                        sheet = "Subject Data") %>% 
  mutate(rawid = toupper(`Subject ID`),
         method = "taqman") %>% 
  select(rawid, 
         nqo1_rs1800566 = `NQO1 genotype: 4559C>T; Pro187Ser; chr16:69711242; rs1800566`,
         method) %>% 
  distinct() %>% 
  group_by(rawid) %>% 
  mutate(row = row_number(),
         method = factor(method),
         nqo1_rs1800566 = factor(nqo1_rs1800566, levels = c("CC", "CT", "TT"), labels = c("0", "1", "2"))) %>% 
  select(-row)

table(taqman$nqo1_rs1800566)

#### merge three datasets (taqman, raw, and imputed genotypes)
vip = raw_clean %>% 
  full_join(imputed) %>% 
  complete(rawid, method) %>% 
  full_join(taqman) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = nqo1_rs1800566) %>% 
  group_by(rawid) %>% 
  fill(gwas_raw, .direction = "downup") %>% 
  fill(gwas_imputed, .direction = "downup") %>% 
  fill(taqman, .direction = "downup") %>% 
  group_by(rawid)  %>% 
  na_if(., "NA") %>% 
  select(-row)%>% 
  distinct()

table(vip$gwas_raw, vip$gwas_imputed, vip$taqman)


write_delim(vip, "results/datasets/tucson_nqo1.tsv", delim = "\t")

rm(imputed)
rm(raw_clean)
rm(taqman)
rm(vip)


############## 
#### PR #####
##############

#### NQO1*2 (rs1800566)
#### GRCh37.p13 chr 16	69745145 G>A
#### GRCh38.p13 chr 16	69711242 G>A

#### rs1800566 in rawdata.bim
#### chr16:69711242:G:A   in phased_autosomes_lifted.bim 

#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/imputed/phased_autosomes_lifted --snp chr16:69711242:G:A  --recode --out nqo1_rs1800566_imputed_sj
#### plink1.9 --bfile ~/../shared/cohort_data/warf_data/gwas/PR/rawdata/rawdata --snp rs1800566  --recode --out nqo1_rs1800566_raw_sj
#### scp steiner@karneslab.pharmacy.arizona.edu:/home/steiner/GWAS_MAY2021/candidate_genes/cyp2c9_2_3_*_tucson.ped ~/Documents/dissertation/GWAS/data/

#### load raw GWAS data 
raw = read_table("results/datasets/nqo1_rs1800566_raw_sj.ped",
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
imputed = read_table("results/datasets/nqo1_rs1800566_imputed_sj.ped",
                          col_names = F) %>% 
  select(id = X2,
         snp1 = X7, 
         snp2 = X8)%>% 
  unite("vip", starts_with("snp")) %>%
  mutate(vip = factor(vip, levels = c("A_A", "A_G", "G_G"), labels = c( "2", "1", "0")), 
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
         vip = `NQO1 genotype: 4559C>T; Pro187Ser; chr16:69711242; rs1800566`
          ) %>% 
  mutate(id = sub(".*?_", "", rawid),
         method = "dmet_luminex",
         vip = factor(vip, levels = c("C/C", "C/T", "T/T"), labels = c("0", "1", "2"))) %>% 
  select(-id)

rm(UPR_ids)
table(UPR$vip)

#### merge three datasets (dmet/luminex, raw, and imputed genotypes)
vip = raw %>% 
  full_join(imputed) %>% 
  complete(rawid, method) %>% 
  full_join(UPR) %>% 
  complete(rawid, method) %>% 
  group_by(rawid) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider( names_from = method, 
               values_from = vip) %>% 
  group_by(rawid) %>% 
  fill(gwas_raw, .direction = "updown") %>% 
  fill(gwas_imputed, .direction = "downup") %>% 
  fill(dmet_luminex, .direction = "downup") %>% 
  group_by(rawid)  %>% 
  na_if(., "NA")%>% 
  select(-row) %>% 
  distinct()

table(vip$gwas_raw, vip$gwas_imputed, vip$dmet_luminex)

write_delim(vip, "results/datasets/sj_nqo1.tsv", delim = "\t")

rm(imputed)
rm(raw)
rm(UPR)
rm(vip)
