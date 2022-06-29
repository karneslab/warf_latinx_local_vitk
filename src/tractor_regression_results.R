#### IN: output of hail.py
#### OUT: plink query
#### as a reminder the model is mt.pheno.dose,[1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.hapcounts1.x, mt.anc1dos.x, mt.anc2dos.x]))
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-03-04
#### Last updated: 2022-06-27



#### load libraries
library(tidyverse)
library(cowplot)


#### load data
lmResults = list.files('results/datasets',
                       pattern = '06-20_vitk90_maf05',
                       full.names = T)
lmResults = lmResults[lmResults != "*.png"]
lmResults


cleanTractor = function(fileName) {
  #### load data
  lm = read_tsv(fileName) %>%
    filter(!is.na(beta))

  clean_hail_tractor <- function(data) {
    data %>%
      dplyr::select(beta, standard_error, t_stat, p_value) %>% ## we need to separate out the values listed in these columns
      names() %>%
      map(function(x) {
        data %>%
          dplyr::select(x) %>%
          separate(x, # separate the columns from beta to pvalue into intercept, hap0, anc0, hap1, anc1, anc2
                   into  = paste0(
                     c("intercept", "hap0", "anc0", "hap1", "anc1", "anc2"),
                     ".",
                     x
                   ),
                   sep = ",")
      }) %>%
      bind_cols() %>%
      bind_cols(data %>% dplyr::select(-beta,-standard_error,-t_stat,-p_value)) %>% ## bind original data (except rows we just separated)
      mutate_at(vars(starts_with("intercept")), ~ gsub("[", "", ., fixed = T)) %>% ## remove all [ signs
      mutate_at(vars(starts_with("anc2")), ~ gsub("]", "", ., fixed = T)) %>% ## remove all ] signs
      separate(variant, into = c("CHR", "BP", "REF", "ALT")) %>%
      mutate_at(vars(-REF,-ALT), as.numeric) %>%
      mutate(match = if_else((anc0.beta > 0 & # this match determines if all ancestries are in the same direction
                                anc1.beta > 0 & ## I never filtered on this but it could be useful 
                                anc2.beta > 0) | (anc0.beta < 0 & anc1.beta < 0 & anc2.beta < 0),
                             "yes",
                             "no"
      ))
  }


  #### clean data
  lm_dat = clean_hail_tractor(lm)


}


tractordats <- lapply(lmResults, cleanTractor)
az = tractordats[[2]]
pr = tractordats[[1]]

tophits_az = az %>% 
  dplyr::select(CHR, BP, REF, ALT, ends_with(c("beta", "p_value")) & starts_with("anc")) %>% 
  arrange(across(ends_with("p_value"))) %>% 
  slice_head(n = 100)

tophits_pr = pr %>% 
  dplyr::select(CHR, BP, REF, ALT, ends_with(c("beta", "p_value")) & starts_with("anc")) %>% 
  arrange(across(ends_with("p_value"))) %>% 
  slice_head(n = 100)

print("Your data is clean :)")



### find replicate signals under 0.0125 in any ancestry 
### the significant ancestry must remain the same in both cohorts

replications = az %>%
  filter(if_any(.cols = c("anc0.p_value", "anc1.p_value", "anc2.p_value"),
                   ~ . < 0.0125)) %>% 
  inner_join(pr  %>%
               filter(if_any(.cols = c("anc0.p_value", "anc1.p_value", "anc2.p_value"),
                        ~ . < 0.0125)), by = c("BP", "CHR", "REF", "ALT"), suffix = c("?tucson", "?sj")) %>% 
  mutate(CHR = paste0("chr",CHR)) %>%
  unite("SNP",CHR:ALT,  sep = ":") %>%
  dplyr::select(-starts_with("n?"),
         -starts_with("match?"),
         -starts_with("f_stat?")) %>% 
  pivot_longer(cols = -SNP,
               names_to = c("stat", "cohort"),
               names_pattern = "(.*)\\?(.*)") %>% 
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  mutate(cohort = factor(cohort, levels = c("tucson", "sj"),
                         labels = c("Tucson", "San Juan"))) %>% 
  group_by(SNP) %>%
  arrange(.by_group = T) %>% 
  group_by(SNP) %>% 
  mutate(anc_sig_match = if_else(all(anc0.p_value< 0.0125)  , "yes_EUR",  ## trying to figure out if the significant pvalue is the same in both cohorts
                                 if_else(all(anc1.p_value< 0.0125) , "yes_NAT", 
                                         if_else(all(anc2.p_value< 0.0125) ,"yes_AFR", "NO_match")))) %>% 
  filter(anc_sig_match != "NO_match")
  

### filter out results with switching betas between cohorts  

replications2 = replications %>% 
  mutate(sig_anc = if_else(anc0.p_value < 0.0125, "EUR", 
                           if_else(anc1.p_value < 0.0125, "NAT",
                                   if_else(anc2.p_value < 0.0125, "AFR", "sig-issue")))) %>% 
  mutate(test_match = ifelse(anc0.p_value < 0.0125 & sign(anc0.beta) == sign(lag(anc0.beta)), "pass",
                             ifelse(anc1.p_value < 0.0125 & sign(anc1.beta) == sign(lag(anc1.beta)), "pass", 
                                    if_else(anc2.p_value < 0.0125 & sign(anc2.beta) == sign(lag(anc2.beta)), "pass", "fail")))) %>% 
  fill(test_match, .direction = "up") %>% 
  filter(test_match != "fail")

length(unique(replications2$SNP))


######## plink query
replications_top = replications2  %>% 
  dplyr::select(SNP, cohort, ends_with(c("beta", "p_value")) & starts_with("anc")) %>% 
  group_by(SNP, cohort) 

replications_out = replications2 %>% 
  dplyr::select(SNP) %>% 
  mutate(SNP = gsub("chr", "", SNP)) 

replications_out

replications_out %>% 
  write_tsv("results/datasets/vitk90_LAadj_maf5_replications.txt", col_names = T)


#### write out summary stats
replications2 %>% 
  write_tsv("results/datasets/vitk90_LAadj_maf5_regression_replicates.tsv")



