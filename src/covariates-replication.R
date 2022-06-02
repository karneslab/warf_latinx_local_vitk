#### UPR Covar file from IWPC dataset
#### complete 1) lai_global_f.py, 2) cyp.R, 3) vkor.R, 
#### 4) localancestry_warf_vips.R, and 5) ancestryqc_absolutepath.sh and pca.R prior to running
#### heidiesteiner@email.arizona.edu
#### 30 Jan 2022 
#### last updated FEB 1 2022



#### load packages 
library(readxl)
library(tidyverse)
library(fastDummies)

#### save vip possible categories at names in list
vips = c('dmet_luminex', 'gwas_imputed', 'gwas_raw')


#### load id dictionary
UPR_ids = read_excel("data/ID List.xlsx") %>%
  mutate(FID =  `Sample_ID_GWAS array`,
         IID = `Patient ID Code`) %>%
  select(-...9, -...10, -...11, -Samples_ID_DMET_Plus_Files, -`Sample_ID_GWAS array`)

#### load clinical data 
UPR <- read_excel("data/UPR.xlsx") %>%
  rename(IID = 'Subject ID') %>%
  mutate(IID = gsub("WPRAPR", "WPRA",IID),
         bsa = sqrt((`Height (cm)`*`Weight (kg)`)/ 3600)
) %>%
  inner_join(UPR_ids, by = "IID") %>% 
  select(IID, gender = Gender, race_reported = `Race (reported)`, 
         age = `Age at Time of Consent`,
         height_cm = `Height (cm)`, weight_kg = `Weight (kg)`, 
         indication = `Primary Indication for Warfarin Treatment`,
         diabetes = `Diabetes`, amiodarone = `Amiodarone (Cordarone)`,
         carbamazepine = `Carbamazepine (Tegretol)`, 
         phenytoin = `Phenytoin (Dilantin)`, rifampicin = `Rifampin or Rifampicin`,
         dose = `Therapeutic Dose of Warfarin`,
         comorbidities = `Comorbidities (Disease States)`,
         target = `Target INR`,
         smoke = `Current Smoker`,bsa) %>% 
  mutate(target = factor(target, levels = c("2-2.5", "2-3", "2.5-3.5"),
                         labels = c("other", "2_3", "25_35")),
         drink = 0)
  


#### load extra clinical data 
UPR_add = read_xlsx("data/Cohort PR GWAS 4 samples add.xlsx") %>% 
  mutate(dose = case_when(`Therapeutic Daily Dose of Warfarin (mg)` == "5 mg QD" ~ 35,
                                                    `Therapeutic Daily Dose of Warfarin (mg)` == "7.5 mg QD; 5 mg Tue" ~ 50,
                                                    `Therapeutic Daily Dose of Warfarin (mg)` == "2.5 mg QD; 5 mg Mon, Wed, Fri" ~ 25),
         #  `Amiodarone (Cordarone)` = if_else(`Amiodarone (Cordarone)` == 1, 1, 0),
         amiodarone = as.factor(`Amiodarone (Cordarone)`),
         carbamazepine = as.factor( `Carbamazepine (Tegretol)`),
         phenytoin = as.factor(`Phenytoin (Dilantin)`),
         rifampicin = as.factor(`Rifampin or Rifampicin`),
         height_cm = 2.54*`Height (inches)`,
         weight_kg = 0.453592*`Weight (lbs)`) %>% 
  select(IID = `Subject ID`, gender = Gender, race_reported = `Self-reported Race`,
         age = `Age at Time of Consent`, height_cm, weight_kg, 
         indication = `Primary Indication for Warfarin Treatment`,
         diabetes = Diabetes, dose, amiodarone, carbamazepine, phenytoin, 
         rifampicin)


#### filter keep ids and clean for downstream 
UPR_keep = UPR_ids %>% 
  select(IID, FID) %>% #### making a list of IDs
  full_join(UPR_add %>% dplyr::select(IID)) %>% #### adding on 4 IDs 
  inner_join(UPR %>% full_join(UPR_add)) #### joining the data to the IDs
  

rm(UPR, UPR_add, UPR_ids)

#### load VKORC1
vkor = read_tsv("results/datasets/covariates/sj_vkor.tsv") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(FID = rawid) %>%
  select(-rawid)


#### merge clinical and vkor
dat_clin_vkor = UPR_keep %>% 
  left_join(vkor) 

visdat::vis_miss(dat_clin_vkor)
rm(vkor, UPR_keep)


#### load CYP2C9 
cyp = read_tsv("results/datasets/covariates/sj_cyp.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(FID = rawid) %>%
  select(-rawid)

#### merge CYP2C9 w/ dat_clin_vkor
dat_clin_vkor_cyp = dat_clin_vkor %>% 
  left_join(cyp)


rm(dat_clin_vkor, cyp)

#### load CYP4F2 
cyp = read_tsv("results/datasets/covariates/sj_cyp4f2.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(FID = rawid) %>%
  select(-rawid) %>% 
  rename_with(.fn = ~paste0("cyp4f2_",.),.cols = any_of(vips))

#### merge CYP2C9 w/ dat_clin_vkor
dat_clin_vkor_cyp_cyp4f2 = dat_clin_vkor_cyp %>% 
  left_join(cyp)


rm(dat_clin_vkor_cyp, cyp)

#### load NQO1
nqo1 = read_tsv("results/datasets/covariates/sj_nqo1.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(FID = rawid) %>%
  select(-rawid) %>% 
  rename_with(.fn = ~paste0("nqo1_",.),.cols = any_of(vips))

#### merge CYP2C9 w/ dat_clin_vkor
dat_clin_vkor_cyp_cyp4f2_nqo1 = dat_clin_vkor_cyp_cyp4f2 %>% 
  left_join(nqo1)


rm(dat_clin_vkor_cyp_cyp4f2, nqo1)

#### load GGCX
ggcx = read_tsv("results/datasets/covariates/sj_ggcx.tsv")%>% 
  select_if(~sum(!is.na(.)) > 0) %>% ## remove completely empty rows
  mutate(FID = rawid) %>%
  select(-rawid) %>% 
  rename_with(.fn = ~paste0("ggcx_",.),.cols = any_of(vips))

#### merge 
dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx = dat_clin_vkor_cyp_cyp4f2_nqo1 %>% 
  left_join(ggcx)


rm(dat_clin_vkor_cyp_cyp4f2_nqo1, ggcx)

#### load local ancestry @ pgx K3
####  Dr. Karnes wants to see vkor_ibs (0,1,2), vkor_nat (0,1,2), vkor_yri (0,1,2) for example
#####
la_vips_K3 = read_csv("results/datasets/covariates/la_warf_vips_K3_pr.csv")%>% 
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

table(la_vips_K3$la_cyp)
table(la_vips_K3$la_cyp_ibs, la_vips_K3$la_cyp_nat, la_vips_K3$la_cyp_yri)

table(la_vips_K3$la_cyp4f2)
table(la_vips_K3$la_cyp4f2_ibs, la_vips_K3$la_cyp4f2_nat, la_vips_K3$la_cyp4f2_yri)

table(la_vips_K3$la_vkor)
table(la_vips_K3$la_vkor_ibs, la_vips_K3$la_vkor_nat, la_vips_K3$la_vkor_yri)

table(la_vips_K3$la_cyp)
table(la_vips_K3$la_cyp_ibs, la_vips_K3$la_cyp_nat, la_vips_K3$la_cyp_yri)

table(la_vips_K3$la_cyp)
table(la_vips_K3$la_cyp_ibs, la_vips_K3$la_cyp_nat, la_vips_K3$la_cyp_yri)

######

#### load local ancestry @ pgx K2
#### need to update the buckets here. Dr. Karnes wants to see vkor_ibs (0,1,2), vkor_nat (0,1,2), vkor_yri (0,1,2) for example
##### 
la_vips_K2 = read_csv("results/datasets/covariates/la_warf_vips_K2_pr.csv")%>% 
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

table(la_vips_K2$la_cyp)
table(la_vips_K2$la_cyp_ibs, la_vips_K2$la_cyp_nat, la_vips_K2$la_cyp_yri)

table(la_vips_K2$la_cyp4f2)
table(la_vips_K2$la_cyp4f2_ibs, la_vips_K2$la_cyp4f2_nat, la_vips_K2$la_cyp4f2_yri)

table(la_vips_K2$la_vkor)
table(la_vips_K2$la_vkor_ibs, la_vips_K2$la_vkor_nat, la_vips_K2$la_vkor_yri)

table(la_vips_K2$la_cyp)
table(la_vips_K2$la_cyp_ibs, la_vips_K2$la_cyp_nat, la_vips_K2$la_cyp_yri)

table(la_vips_K2$la_cyp)
table(la_vips_K2$la_cyp_ibs, la_vips_K2$la_cyp_nat, la_vips_K2$la_cyp_yri)

##### 

#### join la_vips datasets
la_vips = la_vips_K3 %>%
  select(-rawid) %>% 
  full_join(la_vips_K2 %>% select(-rawid)) %>% 
  mutate(FID = IID) %>% 
  select(-IID)

rm(la_vips_K2, la_vips_K3)

#### merge data 
dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la = dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx %>% 
  left_join(la_vips)


rm(la_vips, dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx)

#### load global ancestry 
global = read_table('results/datasets/covariates/sj_lai.txt', col_names = F, skip =1) %>% 
  dplyr::select(plate_id=X1, EUR = X2, AMR = X6, AFR = X10) %>% 
  mutate(plate_id2 = gsub(pattern = "0_0_0_0_", "", x=  plate_id),
         FID = sub(".*?_", "", plate_id2)) ## fix IDs


#### merge 
dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la_ga = dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la %>% 
  left_join(global)


rm(global, dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la)

#### load principal components 
pcas = read_csv("results/datasets/covariates/pr_pcas.csv") %>% 
  rename(FID = IID)

#### merge 
dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la_ga_pca = dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la_ga %>% 
  left_join(pcas)


rm(pcas, dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la_ga)

#### final cleaning 
rep_data = dat_clin_vkor_cyp_cyp4f2_nqo1_ggcx_la_ga_pca %>% 
  mutate(ei = if_else(carbamazepine == 1 | phenytoin == 1 | rifampicin == 1, 1, 0),
         htn = if_else(grepl("hypertension", comorbidities, ignore.case = T),"2", "1"),
         race = if_else(race_reported %in% c("Decline", "Mestizo", "Unknown"), "Mixed/Missing", race_reported),
         race = factor(race, levels = c("black", "Black", "Mixed/Missing", "white", "White"), 
                       labels = c("Black", "Black", "Mixed/Missing", "white", "white")),
         indication = if_else(indication %in% c("4", "4; 1", "4; 3", "6; 4", "3; 4"), "valve", 
                              if_else(indication %in% c("1", "1, 3", "1; 3", "1; 6", "2", "9", "9,8", "9; 3"), "dvt.pe",
                                      if_else(indication %in% c("3", "3; 1", "3; 1; 6", "3; 1; 8", "3; 2", "3; 6", "3; 8"), "af",
                                              if_else(indication %in% c("6", "6; 8"), "tia", "other")))),
         amio = amiodarone,
         ethnicity_latino = "1"
         ) %>% 
  pivot_wider(names_from = "K", values_from = starts_with('la_')) %>%  ## make super wide 
  select(-race_reported, -rifampicin, -carbamazepine, -phenytoin, -comorbidities) %>% 
  dummy_cols(select_columns = c("race", "indication", "target"), ignore_na = T) %>% 
  rename_with( ~ gsub("target_", "target", .x), 
               starts_with("target")) %>% 
  mutate_at(.vars = c("gender", "race_Black", "race_white","race_Mixed/Missing", 
                      "smoke", "drink", "target2_3", "target25_35", "targetother",
                      "indication_dvt.pe",  "indication_tia", "indication_af",
                      "indication_valve", "indication_other", "diabetes", "htn",  "amio", "ei",
                      "ethnicity_latino"), .funs = as_factor) %>% 
  na_if(x = ., y = "NA") %>% 
  mutate_at(vars(contains('la_')), .funs = list(as_factor)) %>% 
  mutate_at(vars(contains(vips)), .funs = list(as_factor)) %>% 
  select(-indication, -amiodarone, -target)
  


visdat::vis_dat(rep_data)


#### write out data 
rep_data %>% 
  select(-plate_id2) %>% 
  write_delim("results/datasets/pr_covariates_31_JAN_22.txt", delim = "\t")

#### write out hail covariate file 
rep_data %>% 
  select(plate_id, dose) %>% 
  write_delim(., "results/datasets/tractor3a_covariates_pr.txt", delim = "\t")


#### write out tractor covariate file https://github.com/Atkinson-Lab/Tractor/wiki/Step-3b:-Tractor-GWAS-(Local)
rep_data %>% 
  select(IID = plate_id, y = dose) %>% 
  write_delim(., "results/datasets/tractor3b_covariates_pr.txt", delim = "\t")

#### write out plink covariate file
rep_data %>% 
  mutate(FID = 0) %>% 
  select(FID, IID = plate_id,  dose, PC1, PC2, PC3) %>% 
  write_delim(., "results/datasets/plink_covariates_pr.txt", delim = "\t")

