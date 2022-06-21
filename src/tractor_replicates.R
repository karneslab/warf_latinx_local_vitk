#### IN: PLINK files of genotypes at replicated sites in tractor_regression_results.R 
#### OUT: FIGURE LA-ADJ REPLICATIONS ASSOCATION W/ DOSE, IWPC adjusted summary stats
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### created 2022-03-06
#### updated 2022-06-20

#### run the below plink code to produce the results for this script
### need to first create new rsids on autosomes files
# for i in {0..2}; do cut -f4,5 --complement phased_autosomes_lifted.80.anc${i}.hapcount.txt > 80.anc${i}.hapcount.hail.txt; done
# plink1.9 --vcf ../../data/imputed/topmed/phased_autosomes_lifted.90.rsid37.vcf.gz --extract vitk90_LAadj_maf5_replications.txt --out vitk90_LAadj_replications --make-bed --const-fid


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
plink = read.plink(bed= "results/datasets/vitk90_LAadj_replications.bed",
                   #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                   #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                   na.strings = c("0", "-9"))


#### convert from plink to dataframe
plink_df = plink[["genotypes"]] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(IID = rowname) %>% 
  rename_with(.cols = c(c(everything(), -IID)), ~ paste0("chr", .x))


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
plink_pr = read.plink(bed= "results/datasets/pr_vitk90_LAadj_replications.bed",
                      #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                      #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                      na.strings = c("0", "-9"))


#### convert from plink to dataframe
plink_df_pr = plink_pr[["genotypes"]] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(IID = rowname)%>% 
  rename_with(.cols = c(c(everything(), -IID)), ~ paste0("chr", .x))


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
reps = results%>% 
  mutate(cohort = "Tucson") %>% 
  rbind(results_pr %>% mutate(cohort = "San Juan"))

reps %>% 
  write_tsv("results/datasets/vitk_tractor_replications_iwpc.tsv")
  



#### plot

#### PLOT NEEDS GENES! 

pheno_geno %>% 
  select(IID, dose, starts_with("chr")) %>% 
  mutate(study = "Tucson, AZ") %>% 
  rbind(pheno_geno_pr  %>% 
          dplyr::select(IID, dose, starts_with("chr")) %>% 
          mutate(study = "San Juan, PR")) %>% 
  mutate(study = factor(study, levels = c("Tucson, AZ", "San Juan, PR"))) %>% 
  pivot_longer(cols = starts_with("chr"), names_to = "SNP") %>% 
  ggplot(aes(x=  value, y = dose, group = SNP, color = SNP)) +
  geom_jitter(alpha =.2, size = .6, width = .1, height = .05) + 
  stat_smooth(geom='line', se=FALSE, na.rm = T, span = 50,
              position = position_jitter(seed = 90, width = .01, height = 1),
              size = 1) + 
  facet_grid(~study) + 
  labs(x = "Variant Allele Copies", 
       y  = "Warfarin Dose (mg/week)", 
       color = "SNP") + 
  guides(color = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(
    legend.position = c(),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size = "120%")
  ) +
  scale_color_viridis_d()

ggsave(plot = last_plot(),
       "results/plots/tractor_replicates_boxplot.png",
       width = 6,
       height = 6,
       unit = "in")
