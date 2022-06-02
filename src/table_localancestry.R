#### Ancestry Results
#### heidi steiner
#### heidiesteiner@email.arizona.edu
#### 2022-03-05


#### load packages 
library(arsenal)
library(tidyverse)



#### load functions
`%!in%` = negate(`%in%`)

#### save vip possible categories at names in list
vips = c( 'cyp', 'cyp4f2', 'vkor', 'nqo1', 'ggcx')

#### load covariate data 
warf_data = read_table("results/datasets/covariates_30JAN22.txt") %>% 
  filter(!grepl("warfer034", IID, ignore.case = T))%>%  # dose outlier removed from analyses
  mutate(study = "tucson",
         gender = factor(gender, levels = c("1", "2"), labels = c("Women", "Men")),
         vkor = factor(vkor_1639_gwas_imputed, levels = c(2, 1,0), labels = c("C/C", "T/C", "T/T")),
         cyp = factor(cyp_stars_gwas_raw, levels = c("*1/*1", "*1/*2", "*1/*3", "*2/*3"), 
                      labels =  c("Extensive Metabolizer", "Intermediate Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer")),
         amio = factor(amio, levels = c("0", "1", "NA"), labels = c("0", "1", "0"))) %>% 
  dplyr::select(starts_with("la_", ignore.case = T)& contains(c("ibs", "yri", "nat")) & contains("_3"), IID, study)

pr_data = read_table("results/datasets/pr_covariates_31JAN22.txt") %>% 
  mutate(study = "pr",
         gender = factor(gender, levels = c("female", "male"), labels = c("Women", "Men")),
         vkor = factor(vkor_1639_gwas_imputed, levels = c("2", "1", "0"), labels = c("T/T", "T/C", "C/C")),
         cyp = factor(cyp_stars_gwas_imputed, levels = c("*1/*1", "*1/*2", "*1/*3", "*2/*3"), 
                      labels =  c("Extensive Metabolizer", "Intermediate Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer")),
         amio = factor(amio, levels = c("0", "1", "N/A"),
                       labels = c("0", "1", "0"))) %>% 
  dplyr::select(starts_with("la_", ignore.case = T) & contains(c("ibs", "yri", "nat")) & contains("_3"), IID, study)


#### merge data frames
t1_data = warf_data %>% 
  full_join(pr_data) %>% 
  pivot_longer(cols = c(-IID, -study), 
               names_to = c("gene", ".value"),
               names_prefix = "la_",
               names_sep = "_") %>% 
  mutate(gene = factor(gene, 
                       levels = c("cyp", "vkor", "cyp4f2", "ggcx", "nqo1"),
                       labels = c("CYP2C9","VKORC1", "CYP4F2", "GGCX", "NQO1")),
         study = factor(study,
                        levels = c("tucson", "pr"),
                        labels = c("Tucson, AZ", "San Juan, PR"))
         ) %>% 
  mutate_at(vars(c(ibs, nat, yri)), 
            list(~ factor(., levels = c("0", "1", "2", NA),
                           labels = c("0", "1", "2")))
            ) %>% 
  dplyr::select(study, gene, IBS = ibs, NAT = nat, YRI = yri) 

t1_data_man = t1_data %>% 
  filter(gene %in% c("CYP2C9", "VKORC1"))
                     
t1_data_sup = t1_data %>%
  anti_join(t1_data_man)

#### create table1

t1 <- tableby(study ~ IBS + NAT + YRI, 
                data=t1_data_man, 
                strata = gene,
                numeric.stats = c("N", "countpct"), 
                test = T)

tab1 = summary(t1, text = TRUE)
tab1

#### print out table 
#### copy and paste this result into Numbers/excel and update to resemble a table

tab1 %>% 
  as.data.frame() %>% 
  rename(Gene = gene, Copies = 2, `Tucson, AZ (n = 141)` = `Tucson, AZ (N=282)`,
         `San Juan, PR (n =  96)` = `San Juan, PR (N=192)`, `Total (n = 237)` = `Total (N=474)`,
         `P value` = `p value`) %>% 
  filter(Copies != "-  N-Miss") %>% 
  mutate(Copies = case_when(Copies == "-  0" ~ "0",
                            Copies == "-  1" ~ "1",
                            Copies == "-  2" ~ "2",
                            Copies == "IBS" ~ "IBS",
                            Copies == "NAT" ~ "NAT",
                            Copies == "YRI" ~ "YRI",
                            TRUE ~ ' ')) %>% 
  write_csv("results/tables/localancestry_cyp_vkor.csv")



#### create table 2 
t2 <- tableby(study ~ IBS + NAT + YRI, 
              data=t1_data_sup, 
              strata = gene,
              numeric.stats = c("N", "countpct"), 
              test = T)

tab2 = summary(t2, text = TRUE)
tab2

#### print out table 
#### copy and paste this result into Numbers/excel and update to resemble a table

tab2 %>% 
  as.data.frame() %>% 
  rename(Gene = gene, Copies = 2, `Tucson, AZ (n = 141)` = `Tucson, AZ (N=423)`,
         `San Juan, PR (n = 96)` = `San Juan, PR (N=288)`, `Total (n = 237)` = `Total (N=711)`,
         `P value` = `p value`) %>% 
  filter(Copies != "-  N-Miss") %>% 
  mutate(Copies = case_when(Copies == "-  0" ~ "0",
                            Copies == "-  1" ~ "1",
                            Copies == "-  2" ~ "2",
                            Copies == "IBS" ~ "IBS",
                            Copies == "NAT" ~ "NAT",
                            Copies == "YRI" ~ "YRI",
                            TRUE ~ ' ')) %>% 
  write_csv("results/tables/localancestry_supp.csv")
