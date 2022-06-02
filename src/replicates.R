#### Explore SNPs that replicated between Tucson and San Juan
#### follows regression_results.R
#### this script plots! 
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Created 2022-03-13

#### run the example below plink code to produce the results for this script
# plink1.9 --bfile ../data/imputed/topmed/phased_autosomes_lifted --extract vitk_pcadj_replications.txt  --out vitk_pcadj_replications --make-bed

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
plink = read.plink(bed= "results/datasets/vitk_pcadj_replications.bed",
                   #bim = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.bim",
                   #  fam = "/Users/heidisteiner/WORK/Warfarin/GWAS/Candidate_Genes/tractorreplications.fam",
                   na.strings = c("0", "-9"))

hwe = fread("results/datasets/vitk_pcadj_replications.hwe") 


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
  mutate(across(where(is.numeric), round, digits = 2))

#### save results
results = models %>% 
  filter(grepl("chr", term))


################ PR ####################


#### load data
plink_pr = read.plink(bed= "results/datasets/pr_vitk_pcadj_replications.bed",
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

hwe_pre = fread("results/datasets/pr_vitk_pcadj_replications.hwe")


#### merge covariates and genotypes
pheno_geno_pr =  pheno_pr %>% 
  select(IID, plate_id,dose, age, bsa, amio, ei, cyp_stars_gwas_imputed, vkor_1639_gwas_imputed) %>% 
  inner_join(plink_df_pr, by = c("plate_id" = "IID")) %>% 
  mutate_if(grepl("chr", colnames(.)), as.numeric) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~.-1)) %>% 
  mutate_if(grepl("chr", colnames(.)), list(~factor(.))) %>% 
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

#### create the dataset with IWPC adjusted summary statistics
#### for ALL the Primary analysis replicated variants
hits = results %>% 
  mutate(cohort = "Tucson",
         SNP = str_sub(term, end = -2)) %>% 
  rbind(results_pr %>% mutate(cohort = "San Juan",
                              SNP = str_sub(term, end = -2))) %>% 
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
                        levels = c(2,1,0),
                        labels = c("0","1","2"))) %>% 
  separate(SNP, into = c("chr", "pos"),remove = F) %>% 
  mutate(SNP = as.factor(SNP), 
         rsid = fct_recode(SNP,
                          rs352892 = "chr2.143000356",
                          rs352889 = "chr2.143002122",
                          rs9688138 = "chr5.81377835",
                          rs10942262 = "chr5.81379677",
                          rs72819791 = "chr6.10493611",
                          rs2445971 = "chr6.69706161" ,
                          `chr6:69718121:A:T` = "chr6.69718121.", 
                          rs2445972 = "chr6.69761455" ,
                          rs72691561 = "chr9.21813242",
                          rs72691562 = "chr9.21813304" ,
                          rs72691563="chr9.21813496",
                          rs72691564 = "chr9.21813519",
                          rs7068224 = "chr10.26271173",
                          rs79928732 = "chr10.69273916",
                          rs28422950 = "chr10.69276516",
                          rs619297 = "chr10.95197958",
                          rs656155 = "chr10.95228249",
                          rs7934133 = "chr11.86559043",
                          rs4943945 = "chr11.86560403",
                          rs4792446 = "chr17.14206724",
                          rs4646342 = "chr17.17589958"),
         gene = fct_recode(rsid,
                           KYNU = "rs352892", 
                           KYNU = "rs352889",
                           ACOT12 = "rs9688138", 
                           ACOT12 = "rs10942262", 
                           `NA` =  "rs72819791", 
                           LMBRD1 =  "rs2445971",
                           LMBRD1 =  "chr6:69718121:A:T",
                           LMBRD1 =  "rs2445972",
                           MTAP =  "rs72691561", 
                           MTAP = "rs72691562",
                           MTAP = "rs72691563",
                           MTAP = "rs72691564",
                           GAD2 =  "rs7068224", 
                           HK1 = "rs79928732",
                           HK1 =  "rs28422950",
                           ACSM6 = "rs619297",
                           LOC107984257  = "rs656155", 
                           ME3 = "rs7934133", 
                           ME3 = "rs4943945", 
                           COX10 = "rs4792446",
                           PEMT = "rs4646342" )) 

#### but you can use the replicates here to decide which SNPs (if any) to highlight? 
h1 = hits %>% 
  group_by(SNP) %>% 
  separate(SNP, into = c("chr", "pos", "a1", "a2"), sep = "\\.", remove = F) %>% 
  mutate(SNP = str_sub(SNP, end = -5) ) %>% 
  inner_join(plot_dat %>% select(gene, chr, pos, rsid), c("chr", "pos")) %>% 
  inner_join(hwe_all %>% 
               mutate(SNP  = gsub(":", ".", SNP),SNP = str_sub(SNP, end = -5)) %>% 
               dplyr::select(SNP, P, cohort = study), by = c("SNP", "cohort")) %>% 
  distinct() %>% 
  unite( "allele", a1:a2) %>% 
  # select(gene, allele, cohort, estimate, std.error, p.value) %>% 
  write_tsv("results/datasets/vitk_pcadj_replications_iwpc.tsv")


  
# plot_dat_labels <- plot_dat %>% 
#   group_by(SNP, value) %>% 
#   mutate(label_y = quantile(dose, .25)) %>% 
#   group_by(SNP) %>% 
#   summarise(
#   label_x = max(as.numeric(value)),
#   label_y = max(as.numeric(label_y))) %>% 
#   distinct()


#### testing this plot because there are more variants to see

plot_dat  %>% 
  ggplot(aes(value, dose, group = SNP, color = gene)) +
  geom_jitter(alpha =.1, size = .4, width = .1, height = .05) +
  geom_smooth( method = "lm", se = F ) + 
  facet_grid(~study) + 
  labs(x = "Variant Allele Copies", 
       y  = "Warfarin Dose (mg/week)", 
       color = "Gene") + 
  guides(color = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(
    legend.position = c()
  )

  
ggsave(plot = last_plot(),
       "results/plots/replicates_lines.png",
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