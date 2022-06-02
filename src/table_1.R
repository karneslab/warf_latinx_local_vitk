#### Table 1 
#### heidiesteiner@email.arizona.edu
#### 30 Jan 2022

#### load packages 
library(tableone)
library(tidyverse)


#### load functions
`%!in%` = negate(`%in%`)


#### load covariate data 
warf_data = read_table("results/datasets/covariates_30JAN22.txt") %>% 
  filter(!grepl("warfer034", IID, ignore.case = T))%>%  # dose outlier removed from analyses
  mutate(study = "tucson",
         gender = factor(gender, levels = c("1", "2"), labels = c("Women", "Men")),
         vkor = factor(vkor_1639_gwas_imputed, levels = c(2, 1,0), labels = c("C/C", "T/C", "T/T")),
         cyp = factor(cyp_stars_gwas_raw, levels = c("*1/*1", "*1/*2", "*1/*3", "*2/*3"), 
                      labels =  c("Extensive Metabolizer", "Intermediate Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer")),
         amio = factor(amio, levels = c("0", "1", "NA"), labels = c("0", "1", "0"))) %>% 
  select(dose, age, bsa, gender, amio, ei, EUR, AMR, AFR, 
         cyp, vkor,study) 

table(warf_data$cyp)

pr_data = read_table("results/datasets/pr_covariates_31_JAN_22.txt") %>% 
  mutate(study = "pr",
         gender = factor(gender, levels = c("female", "male"), labels = c("Women", "Men")),
         vkor = factor(vkor_1639_gwas_imputed, levels = c("2", "1", "0"), labels = c("T/T", "T/C", "C/C")),
         cyp = factor(cyp_stars_gwas_imputed, levels = c("*1/*1", "*1/*2", "*1/*3", "*2/*3"), 
                      labels =  c("Extensive Metabolizer", "Intermediate Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer")),
         amio = factor(amio, levels = c("0", "1", "N/A"),
                       labels = c("0", "1", "0"))) %>% 
  select(dose, age, bsa, gender, amio, ei, EUR, AMR, AFR, 
         cyp, vkor, study)

table(pr_data$cyp)


#### merge data frames
t1_data = warf_data %>% 
  full_join(pr_data) %>% 
  mutate(study = factor(study, levels = c("tucson", "pr"), labels = c("Tucson, AZ", "San Juan, PR")))

#### select varibles for table 1 reporting
all_vars = c("dose", "bsa", "age","gender","EUR", "amio", "ei", "cyp", "vkor", "study")

#### need to table tableone which statistical test to use when parametric test wont work 
fact_vars = c( "gender","cyp", "vkor", "study", "amio", "ei")

#### create table1
t1 = CreateTableOne(vars = all_vars, factorVars = fact_vars, data = t1_data, strata = "study")

#### print out table 
print(t1, quote=F) #### copy and paste this result into Numbers/excel and update to resemble a table
w