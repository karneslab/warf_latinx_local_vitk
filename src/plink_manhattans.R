#####  MANHATTAN PLOTS from assoc files 
##### H Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-02-02

#### example PLINK 
#### plink1.9 --bfile ../phased_autosomes_lifted --pheno /home/steiner/GWAS_MAY2021/covariates.txt --linear --pheno-name sqrtdose --allow-no-sex --remove /home/steiner/GWAS_MAY2021/Tractor/runs/ibs_nat/plink/fail_IDs.txt --maf 0.01 --hwe 1e-6  --out vitk_unadjusted_maf1_hwe6_rmdoseancestryoutliers --ci 0.95 --extract range /home/steiner/GWAS_MAY2021/candidate_genes/vitkgenes_plink.txt --const-fid; done

#### load packages
library(readr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)

#### load data desktop
#### files have the SNP name in GRCh38 and the location is in GRCh37 !!! 


dat_az <- fread("results/datasets/vitk.90.PCadjusted_maf5_hwe6_rmdoseoutliers.assoc.linear") 
dat_sj <- fread("results/datasets/pr_vitk.90.PCadjusted_maf5_hwe6.assoc.linear")


#### load gene reference data
gene_result <- read_delim("data/gene_result.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

#### define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line


don_az <- dat_az %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  filter(CHR <= 23) %>%
  # Add this info to the initial dataset
  left_join(dat_az, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) 
 # add_row(CHR = 1, BP = 40000,BETA =00, P = min(dat_az$P)*.1) ### add a point so the plot y axis is longer

print("CONGRATULATIONS! Data was cleaned")

don_sj <- dat_sj %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  dplyr::select(-chr_len) %>%
  filter(CHR <= 23) %>%
  # Add this info to the initial dataset
  left_join(dat_sj, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) 
 # add_row(CHR = 1, BP = 40000,BETA =00, P = min(dat_az$P)*.1) 

#rm("dat", "gwas_data", "notsig_data", "sig_data")

nans_az = don_az %>% 
  dplyr::select(SNP,BETA, P) %>% 
  filter(is.na(BETA))

nans_sj = don_sj %>% 
  dplyr::select(SNP,BETA, P) %>% 
  filter(is.na(BETA))

don = don_az %>% 
  mutate(cohort = "tucson") %>% 
  full_join(don_sj %>% mutate(cohort = "sj")) %>% 
  mutate(cohort = factor(cohort, 
                         levels = c("tucson", "sj"),
                         labels = c("Tucson, Arizona, USA (n = 141)", "San Juan, Puerto Rico (n = 96)")))

### prepare x axis 
## chromosome name 
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

## make the plot 
p1 = ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=BETA), alpha = .8, size = 1.2) +
  scale_color_gradientn(colors = c("#007FFF", "white", "red"),
                        values = scales::rescale(c(min(don_az$BETA,na.rm = T), -7, 0, 7, max(don_az$BETA, na.rm=T)))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     limits = c(0,6)) +     
  
  facet_wrap(vars(cohort),nrow = 2, scales = "free") + 
  # add genome-wide sig and sugg lines
  # geom_hline(yintercept = -log10(sig), color = "gray40") +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray40") +
  # 
  # Add highlighted points

  # geom_point(data=subset(don, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 15)
  ) +
  labs(x = "Chromosomal Position", y = "-log10(P-value)", color = "Estimate")


p1

ggsave(p1, 
       filename = "results/plots/manhattan_panel_PCadj.png", 
       width = 10, 
       height = 9, 
       units = "in",
       device = "png",
       bg = "transparent"
       )


ggsave(p1, 
       filename = "results/plots/manhattan_panel_PCadj.svg", 
       width = 10, 
       height = 9, 
       units = "in",
       device = "svg",
       bg = "transparent"
)

