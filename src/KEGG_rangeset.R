#### Copy and paste KEGG pathway gene IDs and translate with uniprot
#### Paste uniprot IDs and translate to Ensembl IDs, also keep 
#### This script takes Ensembl IDs to create TSV for filtering genotype data (PLINK, hail)
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### 2022-02-17 


#### load libraries
library(tidyverse)
library(biomaRt)


#### load data 
vitkgene_ids_kegg = read_csv("results/datasets/methods/vitk_kegg_uniprot.csv")
vitkgene_ids_ensemble = read_csv("results/datasets/methods/vitk_uniprot_ensembl.csv")

tot = unique(vitkgene_ids_kegg$`Entry name`)


#### merge data 
vitkgene_ids = vitkgene_ids_kegg %>% 
  rename(From = Entry) %>% 
  inner_join(vitkgene_ids_ensemble) %>% 
  rename(ensembl_gene_id = To)

#### make a list of gene names to pull from biomart
vitkgenes = vitkgene_ids$ensembl_gene_id

#### load biomart 
biolist <- as.data.frame(listMarts())

#### call a specific build here 
ensembl=useMart("ensembl",
                host="https://grch37.ensembl.org"
                )

esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

#### decide which attributes are needed from ensemble 
t2g<-getBM(attributes=c('hgnc_symbol','ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)


vitkgenes <- t2g %>%  
  inner_join(vitkgene_ids) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome_name, start_position, end_position, `Entry name`) %>% 
  mutate(uniprotKB = gsub("_HUMAN", "", `Entry name`)) %>% 
  filter(chromosome_name %in% 1:22) %>% 
  dplyr::select(-`Entry name`)


#### write out PLINK data 
vitkgenes %>% 
  write_tsv("results/datasets/vitkgenes_plink_GRCh37.txt", col_names = F)


#### write out HAIL data 
vitkgenes %>% 
  dplyr::select(contig = chromosome_name, start = start_position, end = end_position) %>% 
  write_tsv("results/datasets/vitkgenes_hail_GRCh37.tsv", col_names = T)
