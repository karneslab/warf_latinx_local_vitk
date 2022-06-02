#### RFMix output --> Local Ancestry at Very Important Warfarin Pharmaocgenes
#### R script to be run in the server/HPC for the Tucson cohort
#### CYP2C9, VKORC1, CYP4F2, NQO1 and GGCX
#### heidiesteiner@email.arizona.edu
#### 29 Jan 2022

#### load libraries
library(tidyverse)
library(plyr)

#### set paths
path_to_bed_files = "ibs_nat_yri/plots/"
name_out_file = "la_warf_vips_K3"

#### read A.bed
ldf <- list()
listbed <- dir(path = path_to_bed_files, pattern="*A.bed")
for (k in 1:length(listbed)){
  ldf[[k]]<-read_tsv(paste(path_to_bed_files,listbed[k], sep = ""),col_names = FALSE)
  ldf[[k]]$iid <- listbed[k]
}

A <- ldply(ldf,rbind) %>% select(-X6,-X8)
colnames(A) <- c("chr","begin","end","gene_pop1","cmbegin","cmend","iid")

A$iid_name<-gsub("0_0_", "", A$iid)
A$iid_name<-gsub("_.*", "", A$iid_name)

#### read B.bed
ldf <- list()
listbed <- dir(path = path_to_bed_files, pattern="*B.bed")
for (k in 1:length(listbed)){
  ldf[[k]]<-read_tsv(paste(path_to_bed_files,listbed[k], sep = ""),col_names = FALSE)
  ldf[[k]]$iid <- listbed[k]
}

B <- ldply(ldf,rbind) %>% select(-X6,-X8)
colnames(B) <- c("chr","begin","end","gene_pop2","cmbegin","cmend","iid")

B$iid_name<-sub("0_0_","",B$iid)
B$iid_name<-sub("_.*","",B$iid_name)

rm(ldf, listbed)

#### CYP2C9
#### 	Chr10: 94938658 - 94990091 GRCh38 
#### 	Chr10: 96698415 - 96749848 GRCh37  

cyp_a <- A %>% 
  filter(chr ==10 )%>%
  filter(begin <= 96698415 & end >= 96749848)%>% 
  select(iid,  CYP2C9_pop1 = gene_pop1) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".A.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


cyp_b <- B %>% 
  filter(chr ==10 )%>%
  filter(begin <= 96698415 & end >= 96749848) %>% 
  select(iid, CYP2C9_pop2 = gene_pop2) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".B.bed", "", iid2, perl = T))%>% 
  select(-iid, -iid2)

#### combine A.bed and B.bed information for gene CYP2C9
cyp <- cyp_a %>% 
  right_join(cyp_b)

table(cyp$CYP2C9_pop1, cyp$CYP2C9_pop2)

rm(cyp_a,cyp_b)

#### VKORC1
####
#### Chr16: 31102175 - 31106118 GRCh37


vkor_a <- A %>% 
  filter(chr ==16 )%>%
  filter(begin <= 31102175 & end >= 31106118)%>% 
  select(iid,  VKORC1_pop1 = gene_pop1) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".A.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


vkor_b <- B %>% 
  filter(chr ==16 )%>%
  filter(begin <= 31102175 & end >= 31106118)%>% 
  select(iid,  VKORC1_pop2 = gene_pop2) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".B.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### combine information from A.bed and B.bed for gene VKORC1

vkor <- inner_join(vkor_a,vkor_b)

table(vkor$VKORC1_pop1, vkor$VKORC1_pop2)

rm(vkor_a, vkor_b)

#### GGCX
#### 	Chr2: 85544720 - 85561527 GRCh38
#### 	Chr2: 85771843 - 85788616 GRCh37


#### scan for gene
ggcx_a <- A %>% 
  filter(chr ==2 )%>%
  filter(begin <= 85771843 & end >= 85788616)%>% 
  select(iid,  GGCX_pop1 = gene_pop1) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".A.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### scan for gene
ggcx_b <- B %>% 
  filter(chr ==2 )%>%
  filter(begin <= 85771843 & end >= 85788616) %>% 
  select(iid,  GGCX_pop2 = gene_pop2) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".B.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)

#### combine A.bed and B.bed information for gene 
ggcx <- inner_join(ggcx_a,ggcx_b)

table(ggcx$GGCX_pop1, ggcx$GGCX_pop2)

rm(ggcx_a, ggcx_b)

#### CYP4F2
#### Chr19 15878023 - 15898074 GRCh38
#### Chr19 15988833 - 16008884 GRCh37 


#### scan for gene
cyp4f2_a <- A %>% 
  filter(chr ==19 )%>%
  filter(begin <= 15988833 & end >= 16008884) %>% 
  select(iid,  CYP4F2_pop1 = gene_pop1) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".A.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### scan for gene
cyp4f2_b <- B %>% 
  filter(chr ==19 )%>%
  filter(begin <= 15988833 & end >= 16008884) %>% 
  select(iid,  CYP4F2_pop2 = gene_pop2) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".B.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### combine A.bed and B.bed information for gene 
cyp4f2 <- inner_join(cyp4f2_a,cyp4f2_b)

table(cyp4f2$CYP4F2_pop1, cyp4f2$CYP4F2_pop2)

rm(cyp4f2_a, cyp4f2_b)

#### NQO1
#### 
#### Chr16 69743304 - 69760463 GRCh37 

#### scan for gene
nqo1_a <- A %>% 
  filter(chr ==16 )%>%
  filter(begin <= 69743304 & end >= 69760463) %>% 
  select(iid,  NQO1_pop1 = gene_pop1) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".A.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### scan for gene
nqo1_b <- B %>% 
  filter(chr ==16 )%>%
  filter(begin <= 69743304 & end >= 69760463) %>% 
  select(iid,  NQO1_pop2 = gene_pop2) %>% 
  mutate(iid2 = sapply(strsplit(iid, split= "_", fixed = TRUE), tail, 1L),
         rawid = gsub(".B.bed", "", iid2, perl = T)) %>% 
  select(-iid, -iid2)


#### combine A.bed and B.bed information for gene 
nqo1 <- inner_join(nqo1_a,nqo1_b)

table(nqo1$NQO1_pop1, nqo1$NQO1_pop2)

rm(nqo1_a, nqo1_b)

#### join all genes
la_vips = cyp %>% 
  left_join(vkor) %>% 
  left_join(nqo1) %>% 
  left_join(ggcx) %>% 
  left_join(cyp4f2)


#### write out data 
write_csv(la_vips, paste(name_out_file,".csv", sep = ""))

