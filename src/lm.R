#### Script to run GWAS in R
#### I'm doing this because PLINK outputs results by increasing dosages and I want by genotype
#### these are PC adjusted !
#### Heidi Steiner
#### heidiesteiner@email.arizona.edu
#### Created 2022-03-11

#### fix file name if rerunning this 

#### load packages 
library(snpStats)
library(tidyverse)
library(broom)

##Pbars...
library(pbapply)
library(pbmcapply)
##

#### load functions
# "%!in%" = negate(`%in%`)

#### local steps
# pheno = list.files('results/datasets', pattern = 'plink_covariates', full.names = T)
# plink_files <- list.files('results/datasets', pattern = 'vitk_maf1_hwe6_rmdoseoutliers.', full.names = T)




# ## find data
plink_files <- list.files('../plink', pattern = 'vitk_maf1_hwe6_rmdoseoutliers.', full.names = T)
plink_files

pheno = list.files('../plink', pattern = 'plink_covariates', full.names = T)
pheno



#### only include PLINK files 
pattern<-c(".bed")

plink <- unique (grep(pattern, plink_files, value=TRUE))
plink

rm(pattern, plink_files)


# fileName = plink[1]
# phenoName = pheno[1]

#### load PLINK data to R 
gwas_lm <- function(fileName, phenoName, multicore = TRUE) {
  dat = read.plink(bed = unique(grep(".bed", fileName, value = T)), na.strings = c("0", "-9"))
  
  #### convert from plink to dataframe
  dat_df = dat[["genotypes"]] %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(IID = rowname)
  
  # dat_df <- dat_df[,c(1, 4500:5000)]
   
  #### load covariates file
  dat_pheno = read_tsv(phenoName) %>%  
    filter(IID != "0_0_WARFER034_1767-JK_Karnes_MEGA_Plate_01_E05_WARFER034") %>%
    as.data.frame()
  
  #### merge covariates and genotypes
  #### this takes awhile could be improved but I think it works
  cat("Finding doses... \n")
  
  #### Find the dose for each ID
  dC <- do.call(rbind, pblapply(dat_df[,1], function(x){
    dat_pheno[which(dat_pheno$IID ==  x), c("dose", "PC1", "PC2", "PC3")]
  }))

  toFactor <- function(x){factor(x, 
                     levels = c("03", "02", "01"),
                     labels = c("00", "01", "02"))
  }
  
  testing <- if(multicore){
    ncores <- detectCores() - 2
    pbmclapply(dat_df[,-c(1)], toFactor,  mc.cores = ncores)
  }else{
    pblapply(dat_df[,-c(1)], toFactor)
  } 
  
  testing <- do.call(cbind.data.frame, testing) ##This can be improved...
  
  pheno_geno <- cbind.data.frame("IID"=dat_df[,1],
                                 dC, 
                                 testing)
  
  cat("Fixing colnames\n")
  
  #### loop regressions
   if(!multicore){ 
    
     colnames(pheno_geno) = gsub(":", ".", colnames(pheno_geno))
     cat("Ready for GWAS\n")
     #### select variables to run regressions on
     myvars = colnames(pheno_geno)
     mySNPs <- myvars[myvars != c("IID", "dose", "PC1", "PC2", "PC3")]
     
     models <- mySNPs %>%
    str_c("dose ~ PC1 + PC2 + PC3 +", .) %>%
    map(.f = ~ lm(formula = as.formula(.x),
                  data = pheno_geno)) %>%
    map(.f = ~ tidy(.x,
                    exponentiate = TRUE,
                    conf.int = TRUE)) %>%
    bind_rows() %>%
    
    #### round all numeric columns
    mutate(across(where(is.numeric), round, digits = 2))
  
  results = models %>%
    filter(grepl("chr", term))
  
  }else{
    ncores <- detectCores() - 2
    targetColumns <- colnames(pheno_geno)[-c(1:2)]
    
    fitRegression <- function(predictor, data=pheno_geno){
      tryCatch({
      predictorFixed = gsub(":", ".", predictor)
      colnames(data)[colnames(data) == predictor] <- predictorFixed
      cat("\r Running GWAS on: ",predictor)
      mod <- summary(lm( formula(paste0("dose ~ PC1+PC2+PC3+", predictorFixed)), data = data))
      tidy(mod,
           exponentiate = TRUE,
           conf.int = TRUE)
      }, error=function(e){
        write(predictor, "err.txt", append = T)
      })
    }
    
    regsRes <-  pbmclapply(targetColumns,
               fitRegression, 
               mc.cores = ncores)
    
    res <- do.call(rbind, regsRes)
    results  = res %>%
      filter(grepl("chr", term))
  }
  
  return(results)
}

#### add this when making the above a function 

print("I keep on running running and running runnning")

summary_stats <- pblapply(plink, gwas_lm, phenoName = pheno)

az_summary_stats = summary_stats[[2]] ### I think these are backwards in server script
pr_summary_stats = summary_stats[[1]]


#### check for replicates
az_summary_stats %>% 
  mutate(cohort = "Tucson") %>% 
  rbind(pr_summary_stats %>% mutate(cohort = "San Juan")) %>% 
  write_tsv("results/datasets/vitk_unadj.tsv")



