#### the output of hail_LAadj.py makes a plink file here
#### to sort for replicates for further analysis
#### as a reminder the model is mt.pheno.dose,[1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.hapcounts1.x, mt.anc1dos.x, mt.anc2dos.x]))
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-03-04
#### Last updated: 2022-03-07


#### load libraries
library(tidyverse)
library(cowplot)

#### load data
lmResults = list.files('results/datasets/',
                       pattern = '2022-03-01',
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
      mutate(match = if_else((anc0.beta > 0 &
                                anc1.beta > 0 &
                                anc2.beta > 0) | (anc0.beta < 0 & anc1.beta < 0 & anc2.beta < 0),
                             "yes",
                             "no"
      )) %>%
      filter(if_any(
        .cols = c("anc0.p_value", "anc1.p_value", "anc2.p_value"),
        ~ . < 0.0125
      )) %>%
      filter(match == "yes")
    
    
  }
  
  
  #### clean data
  lm_dat = clean_hail_tractor(lm)
  
  
}

tractordats <- lapply(lmResults, cleanTractor)
az = tractordats[[1]]
pr = tractordats[[2]]


print("Your data is clean :)")

### find replicate signals

replications = az %>%
  inner_join(pr , by = c("BP", "CHR", "REF", "ALT"),
             suffix = c("?tucson", "?sj")) %>% 
  mutate(CHR = paste0("chr",CHR),
         BP = "17560763") %>% # had to manually change the build info to do the plink query....
  unite("SNP",CHR:ALT,  sep = ":")%>% 
  select(-starts_with("n?"),
         -starts_with("match?"),
         -starts_with("f_stat?")) %>% 
  pivot_longer(cols = -SNP,
               names_to = c("stat", "cohort"),
               names_pattern = "(.*)\\?(.*)") %>% 
  pivot_wider(names_from = stat,
              values_from = value) %>% 
  mutate(cohort = factor(cohort, levels = c("tucson", "sj"),
                         labels = c("Tucson", "San Juan")))

######## plink query
replications  %>% 
  dplyr::select(SNP) %>% 
  write_tsv("results/datasets/tractor_replications.txt", col_names = F)


#### write out summary stats
replications %>% 
  write_tsv("results/datasets/vitk_tractor_regression_replicates.tsv")



