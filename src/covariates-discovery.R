#### Tucson Combining tucson_redcap.tsv with additional covariates
#### complete 1) clean_redcap.R, 2) lai_global_f.py, 3) cyp.R, 4) vkor.R, 
#### 5) localancestry_warf_vips.R, and 6) ancestryqc_absolutepath.sh and pca.R prior to running
#### heidiesteiner@arizona.edu
#### Jan 27 2022
#### last updated 2022-03-10



#### load libraries 
library(tidyverse)
library(visdat)

#### save vip possible categories at names in list
vips = c( 'gwas_imputed', 'taqman', 'gwas_raw')

#### create functions 
`%!in%` = Negate(`%in%`) 

#### load clinical data 
redcap <- read_tsv("results/datasets/covariates/tucson_redcap.tsv")

#### load global ancestry 
global = read_table('results/datasets/covariates/tucson_lai.txt', col_names = F, skip =1) %>% 
  dplyr::select(plate_id=X1, EUR = X2, AMR = X6, AFR = X10) %>% 
  mutate(IID = sapply(strsplit(plate_id, split= "_", fixed = TRUE), tail, 1L)) ## fix IDs

#### merge clinical and global ancestry 
dat_clin_ga = redcap %>% 
  right_join(global) 

rm(redcap, global)

#### load VKORC1 data 
vkor = read_table("results/datasets/covariates/tucson_vkor.tsv") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(IID = rawid) %>%
  select(-rawid)

#### merge clin_ga and VKORC1 genotypes
dat_clin_ga_vkor = dat_clin_ga %>% 
  left_join(vkor)

rm(dat_clin_ga, vkor)

#### load CYP2C9 data 
cyp = read_table("results/datasets/covariates/tucson_cyp.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(IID = rawid) %>% 
  select( -rawid)

#### merge clin_ga_vkor with CYP2C9
dat_clin_ga_vkor_cyp = dat_clin_ga_vkor %>% 
  left_join(cyp) %>% 
  distinct()

rm(dat_clin_ga_vkor, cyp)

#### load CYP4F2 data 
cyp4f2 = read_table("results/datasets/covariates/tucson_cyp4f2.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(IID = rawid) %>% 
  select(-rawid) %>% 
  rename_with(.fn = ~paste0("cyp4f2_",.),.cols = any_of(vips))

#### merge clin_ga_vkor with CYP2C9
dat_clin_ga_vkor_cyp_cyp4f2 = dat_clin_ga_vkor_cyp %>% 
  left_join(cyp4f2) %>% 
  distinct()

rm(dat_clin_ga_vkor_cyp, cyp4f2)

#### load NQO1 data 
nqo1 = read_table("results/datasets/covariates/tucson_nqo1.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(IID = rawid) %>% 
  select( -rawid)

#### merge clin_ga_vkor with CYP2C9
dat_clin_ga_vkor_cyp_cyp4f2_nqo1 = dat_clin_ga_vkor_cyp_cyp4f2 %>% 
  left_join(nqo1) %>% 
  distinct() %>% 
  rename_with(.fn = ~paste0("nqo1_",.),.cols = any_of(vips))

rm(dat_clin_ga_vkor_cyp_cyp4f2, nqo1)

#### load GGCX data 
ggcx = read_table("results/datasets/covariates/tucson_ggcx.tsv")%>% 
  
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(IID = rawid) %>% 
  select(-rawid) %>% 
  rename_with(.fn = ~paste0("ggcx_",.),.cols = any_of(vips))

#### merge clin_ga_vkor with CYP2C9
dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx = dat_clin_ga_vkor_cyp_cyp4f2_nqo1 %>% 
  left_join(ggcx) %>% 
  distinct()

rm(dat_clin_ga_vkor_cyp_cyp4f2_nqo1, ggcx)

#### load principal components 
pcas = read_csv("results/datasets/covariates/warf_pcas.csv")

#### merge clin_ga_vkor 
dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx_pca = dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx %>% 
  left_join(pcas) %>% 
  distinct()

rm(pcas, dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx)

#### load local ancestry @ pgx K3
####  Dr. Karnes wants to see vkor_ibs (0,1,2), vkor_nat (0,1,2), vkor_yri (0,1,2) for example
##### 
la_vips_K3 = read_csv("results/datasets/covariates/la_warf_vips_K3.csv")%>% 
  unite("la_cyp", CYP2C9_pop1, CYP2C9_pop2) %>% 
  unite("la_vkor", VKORC1_pop1, VKORC1_pop2) %>% 
  unite("la_ggcx", GGCX_pop1, GGCX_pop2) %>% 
  unite("la_nqo1", NQO1_pop1, NQO1_pop2) %>% 
  unite("la_cyp4f2", CYP4F2_pop1, CYP4F2_pop2) %>% 
  mutate(IID = rawid, 
         la_cyp_ibs = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                         "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_cyp_nat = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_cyp_yri = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")), 
         la_vkor_ibs = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_vkor_nat = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_vkor_yri = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")), 
         la_cyp4f2_ibs = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_cyp4f2_nat = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_cyp4f2_yri = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")), 
         la_ggcx_ibs = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_ggcx_nat = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_ggcx_yri = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")), 
         la_nqo1_ibs = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_nqo1_nat = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_nqo1_yri = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")),
         K = "3"
         )

##### 

#### load local ancestry @ pgx K2
##### 
la_vips_K2 = read_csv("results/datasets/covariates/la_warf_vips_K2.csv")%>% 
  unite("la_cyp", CYP2C9_pop1, CYP2C9_pop2) %>% 
  unite("la_vkor", VKORC1_pop1, VKORC1_pop2) %>% 
  unite("la_ggcx", GGCX_pop1, GGCX_pop2) %>% 
  unite("la_nqo1", NQO1_pop1, NQO1_pop2) %>% 
  unite("la_cyp4f2", CYP4F2_pop1, CYP4F2_pop2) %>% 
  mutate(IID = rawid, 
         la_cyp_ibs = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("2", "1", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "0", "0",
                                        "1", "0", "0",
                                        "NA")), 
         la_cyp_nat = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "1", "0", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "1", "2", "1",
                                        "0", "1", "0",
                                        "NA")), 
         la_cyp_yri = factor(la_cyp, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                "NA_IBS", "NA_NAT", "NA_IBS",
                                                "IBS_NA", "NAT_NA", "YRI_NA", 
                                                "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                "NA_NA"), 
                             labels = c("0", "0", "1", 
                                        "NA", "NA","NA",
                                        "NA", "NA", "NA",
                                        "0", "0", "1",
                                        "1", "1", "2",
                                        "NA")), 
         la_vkor_ibs = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("2", "1", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "0", "0",
                                         "1", "0", "0",
                                         "NA")), 
         la_vkor_nat = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "1", "0", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "2", "1",
                                         "0", "1", "0",
                                         "NA")), 
         la_vkor_yri = factor(la_vkor, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "0", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "0", "0", "1",
                                         "1", "1", "2",
                                         "NA")), 
         la_cyp4f2_ibs = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                      "NA_IBS", "NA_NAT", "NA_IBS",
                                                      "IBS_NA", "NAT_NA", "YRI_NA", 
                                                      "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                      "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                      "NA_NA"), 
                                labels = c("2", "1", "1", 
                                           "NA", "NA","NA",
                                           "NA", "NA", "NA",
                                           "1", "0", "0",
                                           "1", "0", "0",
                                           "NA")), 
         la_cyp4f2_nat = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                      "NA_IBS", "NA_NAT", "NA_IBS",
                                                      "IBS_NA", "NAT_NA", "YRI_NA", 
                                                      "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                      "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                      "NA_NA"), 
                                labels = c("0", "1", "0", 
                                           "NA", "NA","NA",
                                           "NA", "NA", "NA",
                                           "1", "2", "1",
                                           "0", "1", "0",
                                           "NA")), 
         la_cyp4f2_yri = factor(la_cyp4f2, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                      "NA_IBS", "NA_NAT", "NA_IBS",
                                                      "IBS_NA", "NAT_NA", "YRI_NA", 
                                                      "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                      "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                      "NA_NA"), 
                                labels = c("0", "0", "1", 
                                           "NA", "NA","NA",
                                           "NA", "NA", "NA",
                                           "0", "0", "1",
                                           "1", "1", "2",
                                           "NA")), 
         la_ggcx_ibs = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("2", "1", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "0", "0",
                                         "1", "0", "0",
                                         "NA")), 
         la_ggcx_nat = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "1", "0", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "2", "1",
                                         "0", "1", "0",
                                         "NA")), 
         la_ggcx_yri = factor(la_ggcx, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "0", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "0", "0", "1",
                                         "1", "1", "2",
                                         "NA")), 
         la_nqo1_ibs = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("2", "1", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "0", "0",
                                         "1", "0", "0",
                                         "NA")), 
         la_nqo1_nat = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "1", "0", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "1", "2", "1",
                                         "0", "1", "0",
                                         "NA")), 
         la_nqo1_yri = factor(la_nqo1, levels = c("IBS_IBS", "IBS_NAT","IBS_YRI", 
                                                  "NA_IBS", "NA_NAT", "NA_IBS",
                                                  "IBS_NA", "NAT_NA", "YRI_NA", 
                                                  "NAT_IBS", "NAT_NAT", "NAT_YRI", 
                                                  "YRI_IBS", "YRI_NAT", "YRI_YRI",
                                                  "NA_NA"), 
                              labels = c("0", "0", "1", 
                                         "NA", "NA","NA",
                                         "NA", "NA", "NA",
                                         "0", "0", "1",
                                         "1", "1", "2",
                                         "NA")),
         K = "2"
  )

##### 

#### join la_vips datasets
la_vips = la_vips_K3 %>%
  select(-rawid) %>% 
  full_join(la_vips_K2 %>% select(-rawid))

rm(la_vips_K2, la_vips_K3)

#### merge LA @ Warfarin VIPs with dat_clin_ga_vkor_cyp_pca
dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca = dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx_pca %>% 
  inner_join(la_vips)

rm(dat_clin_ga_vkor_cyp_cyp4f2_nqo1_ggcx_pca, la_vips)

#### clean up for final writing 
dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean = dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca %>% 
  mutate(indication_dvt.pe = if_else(indication_dvt == "1" | indication_pe == "1", 1,0)) %>% 
  mutate_at(.vars = c("sex", "race_AI", "race_white", "smoke", "drink", "target2_3", "target25_35",
                      "targetother", "indication_dvt.pe", "indication_tia", "indicataion_af",
                      "indication_valve", "indication_other", "diabetes", "htn", "met", "amio", "ei",
                      "ethnicity_latino", "vkor_1639_gwas_imputed", "vkor_1639_taqman", "vkor_1173_gwas_raw", 
                      "vkor_1173_gwas_imputed", "cyp_2_gwas_raw", "cyp_2_gwas_imputed", "cyp_3_gwas_raw",
                      "cyp_3_gwas_imputed", "cyp_stars_gwas_raw", "cyp_stars_gwas_imputed","cyp_stars_taqman", "K"
                      ), .funs = as_factor) %>% 
  mutate_at(vars(contains('la_')), .funs = list(as_factor)) %>% 
  mutate_at(vars(contains(vips)), .funs = list(as_factor)) %>% 
  na_if(x = ., y = "NA") %>% 
  droplevels() %>% 
  rename(gender = sex) %>% 
  select(-indication_dvt, -indication_pe) %>% 
  pivot_wider(names_from = "K", values_from = starts_with('la_')) ## make super wide 



#### view the data before outputing anything 
vis_dat(dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean)

rm(dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca)

#### write out clean covarite file :D
write_delim(dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean, "results/datasets/covariates_30JAN22.txt", delim = "\t")
# write_tsv(dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean, "results/datasets/covariates_30JAN22.tsv") #### export to orange w/ tsv

#### write out hail covariate file 
dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean %>% 
  select(plate_id, dose) %>% 
  write_delim(., "results/datasets/tractor3a_covariates_az.txt", delim = "\t")

#### write out tractor covariate file https://github.com/Atkinson-Lab/Tractor/wiki/Step-3b:-Tractor-GWAS-(Local)
dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean%>% 
  select(IID = plate_id, y = dose) %>% 
  write_delim(., "results/datasets/tractor3b_covariates_az.txt", delim = "\t")

#### write out lm covariate file
dat_clin_ga_la_vkor_cyp_cyp4f2_nqo1_ggcx_pca_clean %>% 
  mutate(intercept = 1, FID = 0) %>% 
  select(FID, IID = plate_id, dose, PC1, PC2, PC3) %>% 
  write_delim(., "results/datasets/plink_covariates.txt", delim = "\t")

