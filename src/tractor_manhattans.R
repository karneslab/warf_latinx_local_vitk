#### the output of hail_autosomes_lm.py makes manhattan plots here
#### as a reminder the model is mt.pheno.dose,[1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.hapcounts1.x, mt.anc1dos.x, mt.anc2dos.x]))
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-02-22
#### Last Update: 2022-03-08

#### load libraries
library(tidyverse)
library(cowplot)

#### use these if local else, follow below
# 
# lmResults = list.files('results/datasets', pattern = '2022-03-01', full.names = T)
# 
# gene_result = read_delim("data/gene_result.txt",
#                           "\t", escape_double = FALSE, trim_ws = TRUE)


lmResults <- list.files('results/datasets', pattern = '2022-03-01_vitk_lmunadj', full.names = T)
lmResults = lmResults[lmResults != "*.png"]
lmResults

TractorManhattans <- function(fileName) {
  #### load data
  lm = read_tsv(fileName) %>% ## need to fix the name of this file!
    filter(!is.na(beta))
  
  print("LM results loaded :)")
  gene_result <- read_delim("data/gene_result.txt",
                            "\t",
                            escape_double = FALSE,
                            trim_ws = TRUE)
  
  print("Gene data loaded :)")
  
  #### define significance
  sig = 5e-8 # significant threshold line
  sugg = 1e-6 # suggestive threshold line
  
  
  #### clean data function
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
      bind_cols(data %>% select(-beta,-standard_error,-t_stat,-p_value)) %>% ## bind original data (except rows we just separated)
      # sample_frac(.1) %>%  ## comment out when ready to push to server
      mutate_at(vars(starts_with("intercept")), ~ gsub("[", "", ., fixed = T)) %>% ## remove all [ signs
      mutate_at(vars(starts_with("anc2")), ~ gsub("]", "", ., fixed = T)) %>% ## remove all ] signs
      separate(variant, into = c("CHR", "BP", "REF", "ALT")) %>%
      mutate_at(vars(-REF,-ALT), as.numeric) ## this one doesn't seem to work yet
    
  }
  
  #### clean data
  lm_dat = clean_hail_tractor(lm)
  
  print("Your data is clean :)")
  #### create data for plotting function
  make_ggplot_data <- function(data) {
    data %>%
      
      # Compute chromosome size
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      
      # Calculate cumulative position of each chromosome
      mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
      dplyr::select(-chr_len) %>%
      filter(CHR <= 23) %>%
      # Add this info to the initial dataset
      left_join(data, ., by = c("CHR" = "CHR")) %>%
      
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate(BPcum = BP + tot) %>%
      pivot_longer(
        cols = intercept.beta:anc2.p_value,
        names_to = c(".value", "stat"),
        names_sep = "\\."
      ) %>%
      pivot_longer(cols = starts_with("anc"),
                   names_to = "anc") %>%
      mutate(anc = factor(
        anc,
        levels = c("anc0", "anc1", "anc2"),
        labels = c('IBS', 'NAT', 'YRI')
      )) %>%
      select(-intercept,-hap0,-hap1) %>%
      pivot_wider(names_from = "stat", values_from = "value")
  }
  
  ## create plotting data
  plot_dat = make_ggplot_data(lm_dat)
  
  print("You're ready to plot :)")
  ## stop here - not sure if this is needed. n's sample size and n's variants seem the same per dataset?
  # variants = function(x){
  #   x %>%
  #   group_by(anc) %>%
  #   dplyr::select(starts_with("p_value"), anc) %>%
  #   summarise(variant = sum(!is.na(.)))
  # }
  #
  # variant_n = variants(plot_dat)
  #
  
  ### prepare x axis
  axisdf = plot_dat %>% group_by(CHR) %>% summarize(center = (max(BPcum, na.rm = T) + min(BPcum, na.rm = T)) / 2)
  
  ## make the plot fxn
  anc_manhattans <- function(x) {
    ggplot(x, aes(x = BPcum, y = -log10(p_value))) +
      
      # Show all points
      geom_point(aes(color = beta), alpha = .8, size = 1.2) +
      scale_color_gradientn(colors = c("#007FFF", "thistle2", "red"),
                            values = scales::rescale(c(
                              min(x$beta, na.rm = T),10, 0, 10, max(x$beta, na.rm = T)
                            ))) +
      
      # custom X axis:
      scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
      scale_y_continuous(expand = expansion(mult = c(0, .1)), # no space between xaxis and data but 10% on top
                         breaks = c(0, 2, 4, 6),
                         limits = c(0,6)) +     
      
      # add genome-wide sig and sugg lines
      # geom_hline(yintercept = -log10(sig), color = "gray40") +
      # geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray40") +
      
      # Add highlighted points
      # annotate(geom = "text", x = 2000000000, y = 6.5, label = paste0("n=",variant_n$variant))+
      # geom_point(data=subset(don, is_highlight=="yes"), color="#EF4056", size=2) +
      
      facet_wrap( ~ anc, nrow = 3, 
                  strip.position = "right") +
      
      # Custom the theme:
      theme_minimal() +
      theme(
        legend.position = "none",
        
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(color = "transparent", fill = "transparent"),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
      
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 15),
        
        strip.text = element_text(size = 15),
        strip.background = element_rect(
          color = NA,
          fill = NA,
          size = 1.5,
          linetype = "solid"
          
        )
      ) +
      labs(x = "", y = "", color = "Beta")
  }
  
  plotN <- anc_manhattans(plot_dat)
  
  return(plotN)
}

print("I call the plots :)")
tractormanhattanplots <- lapply(lmResults, TractorManhattans)

cowplot = plot_grid(
  tractormanhattanplots[[1]] ,
  tractormanhattanplots[[2]],
  labels = c('A. Tucson, AZ (n=142)', 'B. San Juan, PR (n=96)'),
  label_x = c(-.12,-.12)
)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  tractormanhattanplots[[1]] +
    theme(rect=element_rect(fill="transparent"),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = col),
          axis.line = element_line(colour = col, color = col, size = 0.1),
          panel.border = element_rect(colour = col),
          panel.background= element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 10))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot_wlegend = plot_grid(cowplot, 
                            legend, 
                            ncol = 1, 
                            rel_heights = c(1, .1))

# cowplot_wlegend

print("Plotting now :)")

ggsave(
  filename = "results/plots/hail_tractor_manhattans_2022-03-01_vitk.png",
  plot = print(cowplot_wlegend),
  height = 9,
  width = 10,
  units = "in",
  device = png,
  bg = "transparent"
)
