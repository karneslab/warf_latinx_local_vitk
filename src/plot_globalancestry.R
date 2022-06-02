#### Plots of RFMIX Global Ancestry Estimate 
#### this comes after 01_covariates.R
#### heidiesteiner@email.arizona.edu
#### 30 Jan 2022

#### ** need to add PR here 

#### load libraries
library(tidyverse)
library(cowplot)

#### load tucson global ancestry 
lai_global = read_table('results/datasets/covariates/tucson_lai.txt', col_names = F, skip =1) %>% 
  dplyr::select(plate_id=X1, EUR = X2, NAT = X6, YRI = X10) %>% 
  mutate(IID = sapply(strsplit(plate_id, split= "_", fixed = TRUE), tail, 1L)) ## fix IDs


#### load covariate file 
pheno = read_tsv("results/datasets/covariates_30JAN22.txt") %>% 
  mutate(IID = gsub("0_", "", IID))


#### create merged data
data = lai_global %>% 
  mutate(ID = str_sub(IID, 9)) %>% 
  inner_join(pheno) 

rm(lai_global, pheno)

#### create long data for plotting 
datalong = data %>% 
  pivot_longer(cols = c( "EUR", "NAT", "YRI"), names_to = "ancestry")

#### create global ancestry stack plot 
p1 = datalong %>%
  ggplot(aes(x=IID, y = value, group = IID)) + 
  
  geom_col(aes(fill = ancestry), width=2)+
  scale_fill_viridis_d(option = "D")+
  
  scale_y_continuous(expand = c(0,0))+
  
  labs(y = "Ancestry Proportion", x = "Individuals")+
  theme(axis.text.x = element_blank( ),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        axis.ticks.x = element_blank())
p1

#### summarize ancestry proportions 
datalong %>% 
  group_by(ancestry) %>% 
  summarise(sum = sum(value)/n())



#### San Juan, Puerto Rico 

#### load gloabl ancestry data 
lai_global_pr <- read_table("results/datasets/covariates/sj_lai.txt", 
                          col_names = FALSE, skip = 1) %>% 
  dplyr::select(IID=X1, EUR = X2, NAT = X6, YRI = X10)

#### load phenotype data 
pheno_pr = read_tsv("results/datasets/pr_covariates_31JAN22.txt") %>% 
  select(-EUR, -AMR, -AFR)

#### create merged data 
data_pr = lai_global_pr %>% 
  inner_join(pheno_pr, by = c("IID" = "plate_id")) 

#### create long data for plotting 
datalong_pr = data_pr %>% 
  pivot_longer(cols = c( "EUR", "NAT", "YRI"), names_to = "ancestry")


#### create stacked global ancestry plot 
p2 = datalong_pr %>%
  ggplot(aes(x=IID, y = value, group = IID)) + 
  geom_col(aes(fill = ancestry), width=2)+
  scale_fill_viridis_d(option ="D")+
  labs( x = "Individuals", fill = "Reference \nPopulations")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank( ),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
p2


#### summarize proportions 
datalong_pr %>% 
  group_by(ancestry) %>% 
  summarise(sum = sum(value)/n())

#### merge cohorts ? 
datalong_merge = datalong %>% 
  mutate(cohort = "Tucson") %>% 
  select(IID, ancestry, value, cohort) %>% 
  full_join(datalong_pr %>% mutate(cohort = "San Juan")) %>% 
  mutate(cohort = factor(cohort, levels = c("Tucson", "San Juan"),
                         labels = c("A. Tucson (n = 142)", "B. San Juan (n = 96)")),
         ancestry = factor(ancestry, levels = c("EUR", "NAT",  "YRI"),
                           labels = c("European", "American", "African")))

#### create a plot that includes all in discovery and replication 
p3 = datalong_merge %>%
  ggplot(aes(x=IID, y = value, group = IID)) + 
  
  geom_col(aes(fill = ancestry), width=2)+
  scale_fill_viridis_d(option ="D")+
  scale_y_continuous(expand = c(0,0))+
  
  labs(y = "Ancestry Proportion", x = "Individuals", fill = "Reference \nPopulations")+
  facet_grid(~cohort, scales = "free")+
  theme(axis.text.x = element_blank( ),
        axis.text.y = element_text(size = 10),
        
        axis.ticks.x = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_text(hjust = -.02, size = 12),
        strip.text.x = element_text(),
        
        axis.title = element_text(size = 15),
        
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15)
        )

p3


#### create plot grid 
cowplot = plot_grid(p1, p2, nrow = 1, rel_widths = 1:1)

#### write out plots 
ggsave( "results/plots/global.png", 
        p3,
        width = 8,
        height = 8,
        unit = "in",
        device = "png")
