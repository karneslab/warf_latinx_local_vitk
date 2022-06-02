#### write out PC1-10 from eigenvec files
#### heidiesteiner@email.arizona.edu
#### Jan 29 2022


#### consider this script doesn't contain PR data 
#### not sure I did PCA on that data...
#### consider CHB population 


#### load libraries

library(tidyverse)

#### Discovery (Tucson, AZ)
#### load data
pcas = read_table("results/datasets/methods/data_snpsonly_clean_chrpos_nodups.hgdp1000ghg19_pos.eigenvec",
                  col_names = F)

#### load 1000genomes data (all pops)
phase3 <- read_table("data/1000GP_Phase3.sample")

table(phase3$POP)


#### select just the reference populations used here
phase3_ref = phase3 %>%  
  filter(POP %in% c("IBS", 
                    "YRI",
                    "CHB")) %>%  #### do I want CHB? 
  arrange(POP)

table(phase3_ref$POP)

rm(phase3)

#### load references
hgdp_1000g <- read_table("results/datasets/methods/hgdp1000ghg19.clean.fam", col_names = F) %>% 
  unite(IID, X1:X2, sep = "_",remove = F) %>% 
  select(IID, X1 = X2)


#### load IDs used in our analysis
hgdp_ref <- read_table("data/hgdp.keep", col_names = F) %>% 
  select(id = X1) %>% 
  mutate(pop = "NAT")

#### join reference data 
ref = phase3_ref %>% 
  select(X1 = ID, POP) %>% 
  full_join(hgdp_ref %>% select(X1 = id, POP=pop)) %>% 
  select(X1, POP) %>% 
  inner_join(hgdp_1000g)

rm(phase3_ref, hgdp_ref, hgdp_1000g)

#### check on NAs
sum(is.na(ref$POP))

table(ref$POP)

#### join PCAs with IDs
dat_long = pcas %>% 
  rename(IID = X2) %>% 
  left_join(ref, by = "IID") %>% 
  select(IID, POP, X3:X12) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "pca") %>% 
  mutate(pca = gsub("X", "", pca),
         pca = as.numeric(pca)-2,
         pca = paste("PC", pca, sep = ""),
         pop = if_else(grepl("warf", IID, ignore.case = T), "study", POP),
         pop = factor(pop, labels = c("Han Chinese (n=103)", "Iberian (n=107)", "Native American (n=46)", "Yoruban (n=108)", "Study (n=142)"), levels = c("CHB", "IBS", "NAT", "YRI","study"))) %>% 
  filter(!is.na(pop)) 

rm(pcas, ref)

dat_wide = dat_long %>% 
  pivot_wider(names_from = "pca", values_from = "value")

#### write out PCs for covariate file
dat_wide %>% 
  filter(grepl("WARF", IID, ignore.case = T)) %>% 
  select(-POP, -pop) %>% 
  mutate(IID = sapply(strsplit(IID, split= "_", fixed = TRUE), tail, 1L)) %>% ## fix IDs
  write_csv("results/datasets/warf_pcas.csv")

rm(dat_wide, dat_long)


#### Replication (San Juan, Puerto Rico)
#### load data
pcas = read_table("results/datasets/methods/pr_data_snpsonly_clean_chrpos_nodups.hgdp1000ghg19_pos.eigenvec",
                  col_names = F)



#### load 1000genomes data (all pops)
phase3 <- read_table("data/1000GP_Phase3.sample")

table(phase3$POP)


#### select just the reference populations used here
phase3_ref = phase3 %>%  
  filter(POP %in% c("IBS", 
                    "YRI",
                    "CHB")) %>%  #### do I want CHB? 
  arrange(POP)

table(phase3_ref$POP)

rm(phase3)

#### load references
hgdp_1000g <- read_table("results/datasets/methods/hgdp1000ghg19.clean.fam", col_names = F) %>% 
  unite(IID, X1:X2, sep = "_",remove = F) %>% 
  select(IID, X1 = X2)


#### load IDs used in our analysis
hgdp_ref <- read_table("data/hgdp.keep", col_names = F) %>% 
  select(id = X1) %>% 
  mutate(pop = "NAT")

#### join reference data 
ref = phase3_ref %>% 
  select(X1 = ID, POP) %>% 
  full_join(hgdp_ref %>% select(X1 = id, POP=pop)) %>% 
  select(X1, POP) %>% 
  inner_join(hgdp_1000g)

rm(phase3_ref, hgdp_ref, hgdp_1000g)

#### check on NAs
sum(is.na(ref$POP))

table(ref$POP)

#### join PCAs with IDs
dat_long = pcas %>% 
  rename(IID = X2) %>% 
  left_join(ref, by = "IID") %>% 
  select(IID, POP, X3:X12) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "pca") %>% 
  mutate(pca = gsub("X", "", pca),
         pca = as.numeric(pca)-2,
         pca = paste("PC", pca, sep = ""),
         pop = if_else(grepl("r*c*", IID, ignore.case = T, perl = T), "study", POP),
         pop = factor(pop, labels = c("Han Chinese (n=103)", "Iberian (n=107)", "Native American (n=46)", "Yoruban (n=108)", "Study (n=142)"), levels = c("CHB", "IBS", "NAT", "YRI","study"))) %>% 
  filter(!is.na(pop)) 

rm(pcas, ref)

#### make wide data for writing 
dat_wide = dat_long %>% 
  pivot_wider(names_from = "pca", values_from = "value")%>% 
  filter(!grepl("HGDP", IID),
         !grepl("_NA", IID),
         !grepl("HG", IID)) %>% 
  select(-POP, -pop) %>% 
  mutate(IID = gsub("0_0_", "", IID),
         IID = sub(".*?_", "", IID))

#### write out PCs for covariate file
dat_wide  %>%
  write_csv("results/datasets/covariates/pr_pcas.csv")


