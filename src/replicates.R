#### IN: PLINK sourced GENOTYPES & HWE of replicated hits from PC-adj models
#### (follows regression_results.R)
#### OUT: Summary stats of replications after IWPC adjustment
#### OUT: FIGURE, linear relationship y=dose,x=variant copies of each SNP in IN files
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Created 2022-03-13
#### Updated 2022-06-19


#### run the EXAMPLE below plink code to produce the results for this script
# plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted.90 --extract vitk_pcadj_replications.txt  --out vitk_pcadj_replications --make-bed

#### Calculate Hardy weinberg at this step as well 
# plink1.9 --bfile vitk_pcadj_replications --hardy --out vitk_pcadj_replications

#### load libraries

library(broom)
library(car)
library(cowplot)
library(data.table)
library(forcats)
library(ggrepel)
library(snpStats)
library(tidyverse)


#### load functions
"%!in%" = negate(`%in%`)

#### load data
plink = read.plink(bed= "results/datasets/vitk90_pcadj_maf05_replications.bed",
                   #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                   #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                   na.strings = c("0", "-9"))

hwe = fread("results/datasets/vitk90_pcadj_maf05_replications.hwe") # spelling error here 


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
  mutate_if(grepl("chr", colnames(.)), list(~factor(.))) %>% 
  mutate_if(grepl("chr", colnames(.)), fct_rev) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
  mutate_at(vars(c(amio, ei, cyp_stars_gwas_raw, vkor_1639_gwas_imputed)), as.factor) 
  

# pheno_geno_hwe = pheno_geno %>% 
#   summarise(across(starts_with("chr"), ~ hwe(as.numeric(.x), "count")$p.lrt)) %>% 
#   pivot_longer()

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
  mutate(across(c(estimate, std.error, statistic), round, digits = 2))

#### save results
results = models %>% 
  filter(grepl("chr", term))


################ PR ####################


#### load data
plink_pr = read.plink(bed= "results/datasets/pr_vitk90_pcadj_maf05_replications.bed",
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

hwe_pre = fread("results/datasets/pr_vitk90_pcadj_maf05_replications.hwe")


#### merge covariates and genotypes
pheno_geno_pr =  pheno_pr %>% 
  select(IID, plate_id,dose, age, bsa, amio, ei, cyp_stars_gwas_imputed, vkor_1639_gwas_imputed) %>% 
  inner_join(plink_df_pr, by = c("plate_id" = "IID")) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~.-1)) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~factor(.))) %>% 
  mutate_if(grepl("chr", colnames(.)), fct_rev) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
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
  mutate(across(c(estimate, std.error, statistic), round, digits = 2))

#### save results
results_pr = models %>% 
  filter(grepl("chr", term)) 

#### create the dataset with IWPC adjusted summary statistics
#### for ALL the Primary analysis replicated variants
hits = results %>% 
  mutate(cohort = "Tucson",
         SNP = str_sub(term, end = -1)) %>% 
  rbind(results_pr %>% mutate(cohort = "San Juan",
                              SNP = str_sub(term, end = -1))) %>% 
  group_by(term) %>% 
  add_count() %>% # finding unmatched SNPs
  filter(n == 2) # remove SNPs that are only represented by one cohort 

length(unique(hits$SNP))


#### plot of the overall SNPs
#### this plot has nothing to do with the IWPC adjustments whatsoever
#### not even sure why it's here 

### first need HWE dataset
hwe_all = hwe %>% 
  mutate(study = "Tucson") %>% 
  full_join(hwe_pre %>% mutate(study = "San Juan"))

plot_dat = pheno_geno %>%  
  mutate(study = "Tucson, AZ") %>% 
  select(IID, dose, which(names(pheno_geno) %in% hits$SNP), study) %>% 
  full_join(pheno_geno_pr %>% mutate(study = "San Juan, PR") %>% 
          select(IID, dose, which(names(pheno_geno) %in% hits$SNP), study)) %>% 
  mutate(study = factor(study, levels = c("Tucson, AZ", "San Juan, PR"))) %>% 
  pivot_longer(cols = starts_with("chr"), names_to = "SNP") %>% 
  mutate(SNP = str_sub(SNP, end = -5),
         value = factor(value,
                        levels = c(1,2,3),
                        labels = c("0","1","2"))) %>% 
  separate(SNP, into = c("chr", "pos"),remove = F) %>% 
  mutate(rsid = as.factor(SNP)) %>% 
  mutate(rsid = fct_recode(rsid, rs619297 = "chr10.95197958",
                           rs650243 = "chr10.95205120",
                           rs8050894 = "chr16.31093188",
                           rs9934438 = "chr16.31093557",
                           rs11078234 = "chr17.14207584"))

# 

#### but you can use the replicates here to decide which SNPs (if any) to highlight? 
hits %>% 
  group_by(SNP) %>% 
  separate(SNP, into = c("chr", "pos", "a1", "a2"), sep = "\\.", remove = F) %>% 
  mutate(SNP = str_sub(SNP, end = -5) ) %>% 
  inner_join(plot_dat %>% select(chr, pos, rsid), c("chr", "pos")) %>% 
  inner_join(hwe_all %>% 
               mutate(SNP  = gsub(":", ".", SNP),SNP = str_sub(SNP, end = -5)) %>% 
               dplyr::select(SNP, P, cohort = study), by = c("SNP", "cohort")) %>% 
  distinct() %>% 
  unite( "allele", a1:a2) %>% 
  # select(gene, allele, cohort, estimate, std.error, p.value) %>% 
  write_tsv("results/datasets/vitk90_maf05_replications_iwpcadj.tsv")


  
# plot_dat_labels <- plot_dat %>% 
#   group_by(SNP, value) %>% 
#   mutate(label_y = quantile(dose, .25)) %>% 
#   group_by(SNP) %>% 
#   summarise(
#   label_x = max(as.numeric(value)),
#   label_y = max(as.numeric(label_y))) %>% 
#   distinct()


#### testing this plot because there are more variants to see

plot_dat   %>% 
  ggplot(aes(x=  value, y = dose, group = SNP, color = rsid)) +
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
       "results/plots/vitk90_PCadj_maf05_replicates_lines.png",
       width = 6,
       height = 6,
       unit = "in")


### archive 


# plot_dat %>% 
#   ggplot(aes(x = value, 
#              y = dose, 
#              group = value)) + 
#   geom_boxplot(aes(fill = value), 
#                outlier.shape = NA,
#                alpha = .5,
#                color = "gray30") +
#   geom_jitter(aes(color = value), 
#               width=0.1,
#               show.legend = F) + 
#   facet_grid(SNP~study) + 
#   labs(x = "Variant Alleles", 
#        y  = "Warfarin Dose (mg/week)", 
#        color = "Variant Copies",
#        fill = "Variant Copies") + 
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank()
#   )