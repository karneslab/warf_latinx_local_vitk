#### the output of lm.R makes manhattan plots on server! 
#### as a reminder the model is dose ~ PC1 + PC2 + PC3 + factor(variant)
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Last Update: 2022-03-15

#### load libraries
library(tidyverse)
library(cowplot)

#### use these if local else, follow below
gene_result = read_delim("data/gene_result.txt",
                          "\t", escape_double = FALSE, trim_ws = TRUE)

#### server data
gene_result <- read_delim("../gene_result.txt", ### fix
                          "\t",
                          escape_double = FALSE,
                          trim_ws = TRUE)


lm3 = read_tsv("results/datasets/vitk_unadj.tsv") %>% 
  mutate(variant = str_sub(term, end = -3)) %>% 
  sample_frac(.2)



#### define significance
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line


#### clean data
dat = lm3 %>% 
  separate(variant, into = c("CHR", "BP", "REF", "ALT"), sep = "\\.") %>% 
  select(CHR, BP, REF, ALT, BETA = estimate, P = p.value, cohort) %>% 
  
  mutate(CHR = gsub("chr","", CHR),
         CHR = as.numeric(as.character(CHR)),
         BP = as.numeric(BP))%>%
  
  arrange(CHR) %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  filter(!grepl("23", CHR)) %>%
  # Add this info to the initial dataset
  left_join(lm3%>% 
              separate(variant, into = c("CHR", "BP", "REF", "ALT"), sep = "\\.") %>% 
              mutate(CHR = gsub("chr","", CHR),
                     CHR = as.numeric(as.character(CHR))), ., by=c("CHR"="CHR")) %>% 
  
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>% 
  mutate_at(c("BP", "tot"), ~as.numeric(.)) %>% 
  mutate( BPcum=BP+tot,
          cohort = factor(cohort, 
                          levels = c("Tucson", "San Juan"),
                          labels = c("Discovery (Tucson, AZ [n = 141])", "Replication (San Juan, PR [n = 96])"))) %>% 
  select(BPcum, P = p.value, CHR, BETA = estimate, term, cohort)

# add_row(CHR = 1, BP = 40000,BETA =00, P = min(dat_az$P)*.1) ### add a point so the plot y axis is longer

print("CONGRATULATIONS! Data was cleaned")


nans = dat %>% 
  dplyr::select(term,BETA, P) %>% 
  filter(!is.na(BETA))

### prepare x axis 
## chromosome name 
axisdf = dat %>% 
  mutate(CHR = as.numeric(as.character(CHR))) %>% 
  arrange(CHR) %>% 
  group_by(CHR) %>% 
  summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

## make the plot 
lm_plot = ggplot(dat, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=BETA), alpha = .8, size = 1.2) +
  scale_color_gradientn(colors = c("#007FFF", "white", "red"),
                        values = scales::rescale(c(min(dat$BETA,na.rm = T), -10, 0, 10, max(dat$BETA, na.rm=T)))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     limits = c(0,6)) +     
  
  
  facet_wrap(vars(cohort),nrow = 2) + 
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 15),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 15)
  ) +
  labs(x = "Chromosomal Position", y = "-log10(P-value)", color = "Estimate")


print("Plotting now :)")
ggsave(
  filename = "results/plots/manhattans_2022-03-15_vitk.png",
  plot = print(lm_plot),
  height = 9,
  width = 10,
  units = "in",
  device = png,
  bg = "transparent"
)
