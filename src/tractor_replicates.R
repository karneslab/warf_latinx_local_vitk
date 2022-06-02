#### Explore SNPs that replicated between Tucson and San Juan in Tractor
#### follows tractor_regression_results.R
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### updated 2022-03-06

#### run the below plink code to produce the results for this script
# plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted 
# --extract tractor_replications.txt --out tractor_vitk_replications 
# --make-bed

#### load libraries
library(broom)
library(car)
library(cowplot)
library(snpStats)
library(tidyverse)


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("snpStats")

#### load functions
"%!in%" = negate(`%in%`)

#### load data
plink = read.plink(bed= "results/datasets/tractor_vitk_replications.bed",
                   #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                   #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                   na.strings = c("0", "-9"))


#### convert from plink to dataframe
plink_df = plink[["genotypes"]] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(IID = rowname)


#### load covariates file 
pheno = read_tsv("results/datasets/covariates_30JAN22.tsv") %>% 
  filter(plate_id != "0_0_WARFER034_1767-JK_Karnes_MEGA_Plate_01_E05_WARFER034") %>% 
  as.data.frame()

#### merge covariates and genotypes 
pheno_geno =  pheno %>% 
  dplyr::select(IID, plate_id, dose, age, bsa, amio, ei, cyp_stars_gwas_raw, vkor_1639_gwas_imputed) %>% 
  inner_join(plink_df, by = c("plate_id" = "IID")) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~.-1)) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~factor(.,
                                                    levels = c("0", "1","2"),
                                                    labels = c("2", "1", "0")))) %>% 
  mutate_if(grepl("chr", colnames(.)), fct_rev) %>% 
  mutate_at(vars(c(amio, ei, cyp_stars_gwas_raw, vkor_1639_gwas_imputed)), as.factor)


  


#### fix column names
colnames(pheno_geno) = gsub(":", ".", colnames(pheno_geno))

#### select variables to run regressions on
myvars = colnames(pheno_geno)[10:ncol(pheno_geno)]

#### loop regressions 
models <- myvars %>%       
  str_c("dose ~ age + bsa + amio + ei +cyp_stars_gwas_raw + vkor_1639_gwas_imputed + ", .) %>%      
  map(                               
    .f = ~lm(                      
      formula = as.formula(.x),     
      data = pheno_geno)) %>%          
  map(
    .f = ~tidy(
      .x, 
      exponentiate = TRUE,       
      conf.int = TRUE)) %>%      
  bind_rows() %>% 
  
  #### round all numeric columns
  mutate(across(where(is.numeric), round, digits = 2))

#### save results
results = models %>% 
  filter(grepl("chr", term))


################ PR ####################


#### load data
plink_pr = read.plink(bed= "results/datasets/pr_tractor_vitk_replications.bed",
                      #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                      #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                      na.strings = c("0", "-9"))


#### convert from plink to dataframe
plink_df_pr = plink_pr[["genotypes"]] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(IID = rowname)


#### covariates
pheno_pr = read_tsv("results/datasets/pr_covariates_31JAN22.txt") %>% 
  as.data.frame()

#### merge covariates and genotypes
pheno_geno_pr =  pheno_pr %>% 
  select(IID, plate_id,dose, age, bsa, amio, ei, cyp_stars_gwas_imputed, vkor_1639_gwas_imputed) %>% 
  inner_join(plink_df_pr, by = c("plate_id" = "IID")) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~.-1)) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~factor(.,
                                                    levels = c("0", "1","2"),
                                                    labels = c("2", "1", "0")))) %>% 
  mutate_if(grepl("chr", colnames(.)), fct_rev) %>% 
  mutate_at(vars(c(amio, ei, cyp_stars_gwas_imputed, vkor_1639_gwas_imputed)), as.factor)

#### fix column names
colnames(pheno_geno_pr) = gsub(":", ".", colnames(pheno_geno_pr))

#### select variables to run regressions on
myvars_pr = colnames(pheno_geno_pr)[10:ncol(pheno_geno_pr)]

#### loop regressions
models <- myvars_pr %>%  
  str_c("dose ~ age+ bsa+amio+ei + cyp_stars_gwas_imputed+vkor_1639_gwas_imputed+", .) %>%    ### edit covarites here
  map(                               
    .f = ~lm(                       # pass the formulas one-by-one to lm()
      formula = as.formula(.x),      # within lm(), the string formula is .x
      data = pheno_geno_pr)) %>%          # dataset
  
  # tidy up each of the glm regression outputs from above
  map(
    .f = ~tidy(
      .x, 
      exponentiate = TRUE,           # exponentiate 
      conf.int = TRUE)) %>%          # return confidence intervals
  
  # collapse the list of regression outputs in to one data frame
  bind_rows() %>% 
  
  # round all numeric columns
  mutate(across(where(is.numeric), round, digits = 2))

#### save results
results_pr = models %>% 
  filter(grepl("chr", term))

#### check for replicates
results%>% 
  mutate(cohort = "Tucson") %>% 
  rbind(results_pr %>% mutate(cohort = "San Juan")) %>% 
  write_tsv("results/datasets/vitk_tractor_replications_iwpc.tsv")
  



#### plot

pheno_geno %>% 
  select(IID, dose, starts_with("chr")) %>% 
  mutate(study = "Tucson, AZ") %>% 
  rbind(pheno_geno_pr  %>% 
          dplyr::select(IID, dose, starts_with("chr")) %>% 
          mutate(study = "San Juan, PR")) %>% 
  mutate(study = factor(study, levels = c("Tucson, AZ", "San Juan, PR"))) %>% 
  pivot_longer(cols = starts_with("chr"), names_to = "SNP") %>% 
  ggplot(aes(x = value, 
             y = dose, 
             group = value)) + 
  geom_boxplot(aes(fill = value), 
               outlier.shape = NA,
               alpha = .5,
               color = "gray30") +
  geom_jitter(aes(color = value), 
              width=0.1,
              show.legend = F) + 
  facet_grid(SNP~study) + 
  labs(x = "Variant Alleles (rs34075240)", 
       y  = "Warfarin Dose (mg/week)", 
       color = "Variant Copies",
       fill = "Variant Copies") + 
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  )



ggsave(plot = last_plot(),
       "results/plots/tractor_replicates_boxplot.png",
       width = 6,
       height = 6,
       unit = "in")
