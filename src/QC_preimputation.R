library(plinkQC)
library(tidyverse)

indir <- "/Users/heidisteiner/WORK/Warfarin/GWAS/warfQC_2021"
name <- 'data'



path2plink = "/anaconda3/bin/plink"
qcdir = "/users/heidisteiner/WORK/Warfarin/GWAS/warfQC_2021"

`%!in%` = negate(`%in%`)
################ per marker
## missingness

fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE)

marker_ids  <- cleanData(
  indir = indir,
  name = name,
  qcdir = indir,
  filterSex = F,
  filterHeterozygosity = F,
  filterSampleMissingness = F,
  filterAncestry = F,
  filterRelated = F,
  filterMAF = F,
  filterHWE = T,
  lmissTh = 0.2,
  hweTh = 0.0000000001,
  mafTh = NULL,
  macTh = NULL,
  path2plink = path2plink,
  verbose = FALSE,
  keep_individuals = NULL,
  remove_individuals = NULL,
  exclude_markers = NULL,
  extract_markers = NULL,
  showPlinkOutput = TRUE
)

name <- 'data_snpsonly.clean'

##### per individual QC 
fail_het_imiss <- check_het_and_miss(indir=indir, qcdir=qcdir, name=name,
                                     interactive=TRUE, path2plink=path2plink,
                                     imissTh = 0.1,
                                     hetTh = 5)


dups = read_table2("data.clean.imiss") %>% 
  mutate(F_MISS = as.numeric(F_MISS)) %>% 
  filter(F_MISS<.2,
         FID %!in% fail_imiss$FID) %>% 
  group_by(IID) %>% 
  add_count() %>% 
  filter(n >1) 

dups2 = dups %>% 
  group_by(IID) %>% 
  slice(which.max(`N_MISS`)) %>% 
  select(FID , IID) 

write_tsv(dups2, "faildups_HS.tsv")

trips = dups %>% 
  filter(n >2,
         FID %!in% dups2$fID) %>% 
  slice(which.max(`N_MISS`)) %>% 
  select(FID,IID)

hapmap = dups %>% 
  filter(grepl("NA", IID)) %>% 
  select(FID,IID)

write_tsv(hapmap, "hapmap_HS.tsv")

write_tsv(trips, "failtrips_HS.tsv")


fail_imiss = fail_het_imiss[["fail_imiss"]] %>% 
  as_tibble() %>% 
  select(FID, IID) %>% 
  as.data.frame()

write_tsv(fail_imiss,"fail_imiss_HS.tsv")
## copy these IDs to fail_IDs.txt
# ggsave(plot = last_plot(), filename = "heterozygosity.png", width = 6)


################################################
##### Relatedness filter function error fix ####
################################################
evaluateDirection <- function(x, y, direction) {
  if (direction == 'ge') x >= y
  else if (direction == 'le') x <= y
  else if (direction == 'gt') x > y
  else if (direction == 'lt') x < y
  else if (direction == 'eq') x == y
  else stop(direction, " as direction in evaluateDirection not known.")
}

makepath <- function(path,file){
  path <- as.list(strsplit(path,'/')[[1]])
  do.call(file.path,c(path,file))
}



relatednessFilter_CRP  = function (relatedness, otherCriterion = NULL, relatednessTh, 
                                   otherCriterionTh = NULL, otherCriterionThDirection = c("gt", 
                                                                                          "ge", "lt", "le", "eq"), relatednessIID1 = "IID1", relatednessIID2 = "IID2", 
                                   relatednessFID1 = NULL, relatednessFID2 = NULL, relatednessRelatedness = "PI_HAT", 
                                   otherCriterionIID = "IID", otherCriterionMeasure = NULL, 
                                   verbose = FALSE) 
{
  if (!(relatednessIID1 %in% names(relatedness))) {
    stop(paste("Column", relatednessIID1, "for relatedness not found!"))
  }
  if (!(relatednessIID2 %in% names(relatedness))) {
    stop(paste("Column", relatednessIID1, "for relatedness not found!"))
  }
  if (!(relatednessRelatedness %in% names(relatedness))) {
    stop(paste("Column", relatednessRelatedness, "for relatedness not found!"))
  }
  iid1_index <- which(colnames(relatedness) == relatednessIID1)
  iid2_index <- which(colnames(relatedness) == relatednessIID2)
  relatedness[, iid1_index] <- as.character(relatedness[, 
                                                        iid1_index])
  relatedness[, iid2_index] <- as.character(relatedness[, 
                                                        iid2_index])
  relatedness_names <- names(relatedness)
  names(relatedness)[iid1_index] <- "IID1"
  names(relatedness)[iid2_index] <- "IID2"
  names(relatedness)[names(relatedness) == relatednessRelatedness] <- "M"
  relatedness_original <- relatedness
  if (!is.null(relatednessFID1) && is.null(relatednessFID2) || 
      is.null(relatednessFID1) && !is.null(relatednessFID2)) {
    stop("Either none or both, relatednessFID1 and relatednessFID2 have to\n             be provided")
  }
  if (!is.null(relatednessFID1) && !is.null(relatednessFID2)) {
    if (!(relatednessFID1 %in% names(relatedness))) {
      stop(paste("Column", relatednessFID1, "for relatedness not found!"))
    }
    if (!(relatednessFID2 %in% names(relatedness))) {
      stop(paste("Column", relatednessFID2, "for relatedness not found!"))
    }
    fid1_index <- which(colnames(relatedness) == relatednessFID1)
    fid2_index <- which(colnames(relatedness) == relatednessFID2)
    names(relatedness)[fid1_index] <- "FID1"
    names(relatedness)[fid2_index] <- "FID2"
    relatedness$FID1 <- as.character(relatedness$FID1)
    relatedness$FID2 <- as.character(relatedness$FID2)
  }
  relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, 
                                ~M)
  sortedIDs <- data.frame(t(apply(relatedness, 1, function(pair) {
    c(sort(c(pair[1], pair[2])))
  })), stringsAsFactors = FALSE)
  keepIndex <- which(!duplicated(sortedIDs))
  relatedness_original <- relatedness_original[keepIndex, 
  ]
  relatedness <- relatedness[keepIndex, ]
  highRelated <- dplyr::filter_(relatedness, ~M > relatednessTh)
  if (nrow(highRelated) == 0) {
    return(list(relatednessFails = NULL, failIDs = NULL))
  }
  uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
  failIDs_other <- NULL
  if (!is.null(otherCriterion)) {
    #  otherCriterionThDirection <- match.arg(otherCriterionThDirection)
    if (!(otherCriterionMeasure %in% names(otherCriterion))) {
      stop(paste("Column", otherCriterionMeasure, "for otherCriterion not found!"))
    }
    if (!(otherCriterionIID %in% names(otherCriterion))) {
      stop(paste("Column", otherCriterionIID, "for otherCriterion not found!"))
    }
    names(otherCriterion)[names(otherCriterion) == otherCriterionMeasure] <- "M"
    names(otherCriterion)[names(otherCriterion) == otherCriterionIID] <- "IID"
    if (any(!uniqueIIDs %in% otherCriterion$IID)) {
      stop("Not all IIDs provided in relatedness are contained in", 
           "otherCriterion")
    }
    fail_other <- apply(highRelated, 1, function(pair) {

      failID1 <- NULL
      failID2 <- NULL
      one <- evaluateDirection(otherCriterion$M[otherCriterion$IID == 
                                                  pair[1]], otherCriterionTh, direction = otherCriterionThDirection)
      two <- evaluateDirection(otherCriterion$M[otherCriterion$IID == 
                                                  pair[2]], otherCriterionTh, direction = otherCriterionThDirection)
      if (one) 
        failID1 <- pair[1]
      if (two) 
        failID2 <- pair[2]
      return(c(failID1, failID2))
    })
    
    failIDs_other <- unique(unlist(fail_other))
    highRelated <- highRelated[!(highRelated$IID1 %in% failIDs_other | 
                                   highRelated$IID2 %in% failIDs_other), ]
    if (nrow(highRelated) == 0) {
      if (verbose) {
        message("Relatedness cannot be evaluated as all individuals ", 
                "involved fail due to otherCriterion")
      }
      return(list(relatednessFails = NULL, failIDs = NULL))
    }
    uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
  }
  
  allRelated <- c(highRelated$IID1, highRelated$IID2)
  multipleRelative <- unique(allRelated[duplicated(allRelated)])
  singleRelative <- uniqueIIDs[!uniqueIIDs %in% multipleRelative]
  highRelatedMultiple <- highRelated[highRelated$IID1 %in% 
                                       multipleRelative | highRelated$IID2 %in% multipleRelative, 
  ]
  highRelatedSingle <- highRelated[highRelated$IID1 %in% singleRelative & 
                                     highRelated$IID2 %in% singleRelative, ]
  if (is.null(nrow(singleRelative)) == F) {
    if (!is.null(otherCriterion)) {
      failIDs_single <- apply(highRelatedSingle, 1, function(pair) {
        print(highRelatedSingle)
        one_two <- evaluateDirection(otherCriterion$M[otherCriterion$IID == 
                                                        pair[1]], otherCriterion$M[otherCriterion$IID == 
                                                                                     pair[2]], direction = otherCriterionThDirection)
        two_one <- evaluateDirection(otherCriterion$M[otherCriterion$IID == 
                                                        pair[2]], otherCriterion$M[otherCriterion$IID == 
                                                                                     pair[1]], direction = otherCriterionThDirection)
        if (one_two) 
          failID <- pair[1]
        else if (two_one) 
          failID <- pair[2]
        else failID <- pair[1]
        return(failID)
      })
    } else {
      failIDs_single <- highRelatedSingle[, 1]
    }
  } else {
    failIDs_single <- NULL
  }
  if (length(multipleRelative) != 0) {
    relatedPerID <- lapply(multipleRelative, function(x) {
      print(multipleRelative)
      tmp <- highRelatedMultiple[rowSums(cbind(highRelatedMultiple$IID1 %in% 
                                                 x, highRelatedMultiple$IID2 %in% x)) != 0, 1:2]
      rel <- unique(unlist(tmp))
      return(rel)
    })
    names(relatedPerID) <- multipleRelative
    relatedPerID = relatedPerID[which(unlist(lapply(relatedPerID, length))>2)]
    #  relatedPerID = relatedPerID[c(1)]
    
    keepIDs_multiple <- lapply(relatedPerID, function(x) {
      print(relatedPerID)
      pairwise <- t(combn(x, 2))
      index <- (highRelatedMultiple$IID1 %in% pairwise[, 
                                                       1] & highRelatedMultiple$IID2 %in% pairwise[, 
                                                                                                   2]) | (highRelatedMultiple$IID1 %in% pairwise[, 
                                                                                                                                                 2] & highRelatedMultiple$IID2 %in% pairwise[, 
                                                                                                                                                                                             1])
      combination <- highRelatedMultiple[index, ]
      combination_graph <- igraph::graph_from_data_frame(combination, 
                                                         directed = FALSE)
      all_iv_set <- igraph::ivs(combination_graph)
      length_iv_set <- sapply(all_iv_set, function(x) length(x))
      print(length_iv_set)
      if (all(length_iv_set == 1)) {
        occurrence <- sapply(x, function(id) {
          print(x)
          sum(sapply(relatedPerID, function(idlist) id %in% 
                       idlist))
        })
        if (length(unique(occurrence)) == 1) {
          nonRelated <- sort(x)[1]
        } else {
          nonRelated <- names(occurrence)[which.min(occurrence)]
        }
      }else {
        nonRelated <- all_iv_set[which.max(length_iv_set)]
      }
      return(nonRelated)
    })
    keepIDs_multiple <- unique(unlist(keepIDs_multiple))
    failIDs_multiple <- c(multipleRelative[!multipleRelative %in% 
                                             keepIDs_multiple])
  } else {
    failIDs_multiple <- NULL
  }
  allFailIIDs <- c(failIDs_single, failIDs_multiple, failIDs_other)
  relatednessFails <- lapply(allFailIIDs, function(id) {
    print(allFailIIDs)
    fail_inorder <- relatedness_original$IID1 == id & relatedness_original$M > 
      relatednessTh
    fail_inreverse <- relatedness_original$IID2 == id & 
      relatedness_original$M > relatednessTh
    if (any(fail_inreverse)) {
      inreverse <- relatedness_original[fail_inreverse, 
      ]
      if (!is.null(relatednessFID2)) {
        id1 <- c(iid1_index, fid1_index)
        id2 <- c(iid2_index, fid2_index)
      } else {
        id1 <- iid1_index
        id2 <- iid2_index
      }
      inreverse[, c(id1, id2)] <- inreverse[, c(id2, id1)]
      names(inreverse) <- relatedness_names
    }else {
      inreverse <- NULL
    }
    inorder <- relatedness_original[fail_inorder, ]
    names(inorder) <- relatedness_names
    return(rbind(inorder, inreverse))
  })
  relatednessFails <- do.call(rbind, relatednessFails)
  if (nrow(relatednessFails) == 0) {
    relatednessFails <- NULL
    failIDs <- NULL
  } else {
    names(relatednessFails) <- relatedness_names
    rownames(relatednessFails) <- 1:nrow(relatednessFails)
    uniqueFails <- relatednessFails[!duplicated(relatednessFails[, 
                                                                 iid1_index]), ]
    if (!is.null(relatednessFID2)) {
      failIDs <- data.frame(FID = uniqueFails[, fid1_index], 
                            IID = uniqueFails[, iid1_index], stringsAsFactors = FALSE)
    } else {
      failIDs <- data.frame(IID = uniqueFails[, iid1_index], 
                            stringsAsFactors = FALSE)
    }
  }
  return(list(relatednessFails = relatednessFails, failIDs = failIDs))
}

evaluate_check_relatedness_CRP <- function(qcdir, name, highIBDTh=0.1875,
                                       imissTh=0.03, interactive=FALSE,
                                       legend_text_size = 5,
                                       legend_title_size = 7,
                                       axis_text_size = 5,
                                       axis_title_size = 7,
                                       title_size = 9,
                                       verbose=FALSE) {
  
  prefix <- makepath(qcdir, name)
  
  if (!file.exists(paste(prefix, ".imiss", sep=""))){
    stop("plink --missing output file: ", prefix,
         ".imiss does not exist.")
  }
  if (!file.exists(paste(prefix, ".genome",sep=""))){
    stop("plink --genome output file: ", prefix,
         ".genome does not exist.")
  }
  testNumerics(numbers=highIBDTh, positives=highIBDTh, proportions=highIBDTh)
  names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
  imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                      as.is=TRUE, stringsAsFactors=FALSE)
  if (!all(names_imiss == names(imiss))) {
    stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
  }
  names_genome <- c("FID1", "IID1", "FID2", "IID2", "RT", "EZ", "Z0", "Z1",
                    "Z2", "PI_HAT", "PHE", "DST", "PPC", "RATIO")
  genome <- read.table(paste(prefix, ".genome", sep=""), header=TRUE,
                       as.is=TRUE, stringsAsFactors=FALSE)
  if (!all(names_genome == names(genome))) {
    stop("Header of ", prefix, ".genome is not correct. Was your
             file generated with plink --genome?")
  }
   fail_highIBD <- relatednessFilter_CRP(relatedness=genome,  otherCriterion=imiss,
                                     relatednessTh=highIBDTh,
                                     relatednessFID1="FID1",
                                     relatednessFID2="FID2",
                                     otherCriterionTh=imissTh,
                                     otherCriterionThDirection="gt",
                                     otherCriterionMeasure="F_MISS" )

  genome$PI_HAT_bin <- ifelse(genome$PI_HAT > 0.05, 0, 1)
   p_allPI_HAT <- ggplot(genome, aes_string('PI_HAT'))
   p_allPI_HAT <- p_allPI_HAT + geom_histogram(binwidth = 0.005,
                                               fill="#66a61e") +
     ylab("Number of pairs") +
     xlab("Estimated pairwise IBD (PI_HAT)") +
     ggtitle("IBD for all sample pairs") +
     geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
     theme_bw() +
     theme(legend.text = element_text(size = legend_text_size),
           legend.title = element_text(size = legend_title_size),
           title = element_text(size = legend_text_size),
           axis.text = element_text(size = axis_text_size),
           axis.title = element_text(size = axis_title_size))
   p_highPI_HAT <- ggplot(dplyr::filter(genome, .data$PI_HAT_bin == 0),
                          aes_string('PI_HAT'))
   p_highPI_HAT <- p_highPI_HAT + geom_histogram(binwidth = 0.005,
                                                 fill="#e6ab02") +
     ylab("Number of pairs") +
     xlab("Estimated pairwise IBD (PI_HAT)") +
     ggtitle("IBD for sample pairs with PI_HAT >0.1") +
     geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
     theme_bw() +
     theme(legend.text = element_text(size = legend_text_size),
           legend.title = element_text(size = legend_title_size),
           title = element_text(size = legend_text_size),
           axis.text = element_text(size = axis_text_size),
           axis.title = element_text(size = axis_title_size))
   p_histo <- cowplot::plot_grid(p_allPI_HAT, p_highPI_HAT)
   title <- cowplot::ggdraw() +
     cowplot::draw_label("Relatedness estimated as pairwise IBD (PI_HAT)",
                         size=title_size)
   p_IBD <- cowplot::plot_grid(title, p_histo, ncol = 1,
                               rel_heights = c(0.1, 1))
   if (interactive) print(p_IBD)
  return(list(fail_highIBD=fail_highIBD$relatednessFails
              ,failIDs=fail_highIBD$failIDs
            , p_IBD=p_IBD,
              plot_data = genome
               ))
 }

check_relatedness_CRP = function (indir, name, qcdir = indir, highIBDTh = 0.1875, genomebuild = "hg19", 
          imissTh = 0.03, run.check_relatedness = TRUE, interactive = FALSE, 
          verbose = FALSE, mafThRelatedness = 0.1, path2plink = NULL, 
          keep_individuals = NULL, remove_individuals = NULL, exclude_markers = NULL, 
          extract_markers = NULL, legend_text_size = 5, legend_title_size = 7, 
          axis_text_size = 5, axis_title_size = 7, title_size = 9, 
          showPlinkOutput = TRUE) {
  if (run.check_relatedness) {
    run <- run_check_relatedness(indir = indir, qcdir = qcdir, 
                                 name = name, verbose = verbose, mafThRelatedness = mafThRelatedness, 
                                 path2plink = path2plink, highIBDTh = highIBDTh, 
                                 keep_individuals = keep_individuals, remove_individuals = remove_individuals, 
                                 exclude_markers = exclude_markers, extract_markers = extract_markers, 
                                 showPlinkOutput = showPlinkOutput)
  }
  fail <- evaluate_check_relatedness_CRP(qcdir = qcdir, name = name, 
                                     highIBDTh = highIBDTh, imissTh = imissTh, interactive = interactive, 
                                     legend_text_size = legend_text_size, legend_title_size = legend_title_size, 
                                     axis_text_size = axis_text_size, axis_title_size = axis_title_size, 
                                     title_size = title_size, verbose = verbose)
  return(fail)
}

perIndividualQC_CRP <- function(indir, name, qcdir=indir,
                                dont.check_sex=FALSE,
                                do.run_check_sex=TRUE, do.evaluate_check_sex=TRUE,
                                maleTh=0.8, femaleTh=0.2,
                                externalSex=NULL, externalMale="M",
                                externalSexSex="Sex", externalSexID="IID",
                                externalFemale="F", fixMixup=FALSE,
                                dont.check_het_and_miss=FALSE,
                                do.run_check_het_and_miss=TRUE,
                                do.evaluate_check_het_and_miss=TRUE,
                                imissTh=0.03, hetTh=3,
                                dont.check_relatedness=FALSE,
                                do.run_check_relatedness=TRUE,
                                do.evaluate_check_relatedness=TRUE,
                                highIBDTh=0.1875,
                                mafThRelatedness=0.1,
                                genomebuild='hg19',
                                dont.check_ancestry=FALSE,
                                do.run_check_ancestry=TRUE,
                                do.evaluate_check_ancestry=TRUE,
                                prefixMergedDataset, europeanTh=1.5,
                                defaultRefSamples = c("HapMap", "1000Genomes"),
                                refSamples=NULL, refColors=NULL,
                                refSamplesFile=NULL, refColorsFile=NULL,
                                refSamplesIID="IID", refSamplesPop="Pop",
                                refColorsColor="Color", refColorsPop="Pop",
                                studyColor="#2c7bb6", label_fail=TRUE,
                                highlight_samples = NULL,
                                highlight_type =
                                  c("text", "label", "color", "shape"),
                                highlight_text_size = 3,
                                highlight_color = "#c51b8a",
                                highlight_shape = 17,
                                highlight_legend = FALSE,
                                interactive=FALSE, verbose=TRUE,
                                keep_individuals=NULL,
                                remove_individuals=NULL,
                                exclude_markers=NULL,
                                extract_markers=NULL,
                                legend_text_size = 5,
                                legend_title_size = 7,
                                axis_text_size = 5,
                                axis_title_size = 7,
                                subplot_label_size = 9,
                                title_size = 9,
                                path2plink=NULL, showPlinkOutput=TRUE) {
  
  missing_genotype <- NULL
  highIBD <- NULL
  outlying_heterozygosity <- NULL
  mismatched_sex <- NULL
  ancestry <- NULL
  
  p_sexcheck <- NULL
  p_het_imiss <- NULL
  p_relatedness <- NULL
  p_ancestry <- NULL
  
  out <- makepath(qcdir, name)
  
  if (!dont.check_sex) {
    if (do.run_check_sex) {
      run <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
                           path2plink=path2plink,
                           showPlinkOutput=showPlinkOutput,
                           keep_individuals=keep_individuals,
                           remove_individuals=remove_individuals,
                           exclude_markers=exclude_markers,
                           extract_markers=extract_markers,
                           verbose=verbose)
    }
    if (do.evaluate_check_sex) {
      if (verbose) {
        message("Identification of individuals with discordant sex ",
                "information")
      }
      fail_sex <- evaluate_check_sex(qcdir=qcdir, indir=indir, name=name,
                                     maleTh=maleTh, femaleTh=femaleTh,
                                     externalSex=externalSex,
                                     externalMale=externalMale,
                                     externalFemale=externalFemale,
                                     externalSexSex=externalSexSex,
                                     externalSexID=externalSexID,
                                     verbose=verbose,
                                     path2plink=path2plink,
                                     showPlinkOutput=showPlinkOutput,
                                     fixMixup=fixMixup,
                                     label_fail=label_fail,
                                     highlight_samples =
                                       highlight_samples,
                                     highlight_type = highlight_type,
                                     highlight_text_size =
                                       highlight_text_size,
                                     highlight_color = highlight_color,
                                     highlight_shape = highlight_shape,
                                     highlight_legend = highlight_legend,
                                     legend_text_size = legend_text_size,
                                     legend_title_size =
                                       legend_title_size,
                                     axis_text_size = axis_text_size,
                                     axis_title_size = axis_title_size,
                                     title_size = title_size,
                                     interactive=FALSE)
      write.table(fail_sex$fail_sex[,1:2],
                  file=paste(out, ".fail-sexcheck.IDs",
                             sep=""),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      if (!is.null(fail_sex$fail_sex) &&
          nrow(fail_sex$fail_sex) != 0) {
        mismatched_sex<- select(fail_sex$fail_sex,
                                .data$FID, .data$IID)
      }
      if (!is.null(fail_sex$mixup)) {
        write.table(fail_sex$mixup[,1:2],
                    file=paste(out, ".sexcheck_mixup.IDs",
                               sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
      }
      p_sexcheck <- fail_sex$p_sexcheck
    }
  }
  if (!dont.check_het_and_miss) {
    if (do.run_check_het_and_miss) {
      run_miss <- run_check_missingness(qcdir=qcdir, indir=indir,
                                        name=name,
                                        path2plink=path2plink,
                                        showPlinkOutput=showPlinkOutput,
                                        keep_individuals=keep_individuals,
                                        remove_individuals=remove_individuals,
                                        exclude_markers=exclude_markers,
                                        extract_markers=extract_markers,
                                        verbose=verbose)
      run_het <- run_check_heterozygosity(qcdir=qcdir, indir=indir,
                                          name=name,
                                          path2plink=path2plink,
                                          showPlinkOutput=showPlinkOutput,
                                          keep_individuals=keep_individuals,
                                          remove_individuals=remove_individuals,
                                          exclude_markers=exclude_markers,
                                          extract_markers=extract_markers,
                                          verbose=verbose)
    }
    if (do.evaluate_check_het_and_miss) {
      if (verbose) {
        message("Identification of individuals with outlying missing ",
                "genotype or heterozygosity rates")
      }
      fail_het_imiss <-
        evaluate_check_het_and_miss(qcdir=qcdir, name=name,
                                    imissTh=imissTh,
                                    hetTh=hetTh, label_fail=label_fail,
                                    highlight_samples =
                                      highlight_samples,
                                    highlight_type = highlight_type,
                                    highlight_text_size =
                                      highlight_text_size,
                                    highlight_color = highlight_color,
                                    highlight_shape = highlight_shape,
                                    highlight_legend = highlight_legend,
                                    legend_text_size = legend_text_size,
                                    legend_title_size =
                                      legend_title_size,
                                    axis_text_size = axis_text_size,
                                    axis_title_size = axis_title_size,
                                    title_size = title_size,
                                    interactive=FALSE)
      write.table(fail_het_imiss$fail_imiss[,1:2],
                  file=paste(out, ".fail-imiss.IDs",
                             sep=""),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      if (!is.null(fail_het_imiss$fail_imiss) &&
          nrow(fail_het_imiss$fail_imiss) != 0) {
        missing_genotype <- select(fail_het_imiss$fail_imiss,
                                   .data$FID, .data$IID)
      }
      
      write.table(fail_het_imiss$fail_het[,1:2],
                  file=paste(out, ".fail-het.IDs",
                             sep=""),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      if (!is.null(fail_het_imiss$fail_het) &&
          nrow(fail_het_imiss$fail_het) != 0) {
        outlying_heterozygosity <- select(fail_het_imiss$fail_het,
                                          .data$FID, .data$IID)
      } else {
        outlying_heterozygosity <- NULL
      }
      p_het_imiss <- fail_het_imiss$p_het_imiss
    }
  }
  if (!dont.check_relatedness) {
    if (do.run_check_relatedness) {
      run <- run_check_relatedness(qcdir=qcdir, indir=indir, name=name,
                                   path2plink=path2plink,
                                   mafThRelatedness=mafThRelatedness,
                                   genomebuild=genomebuild,
                                   showPlinkOutput=showPlinkOutput,
                                   keep_individuals=keep_individuals,
                                   remove_individuals=remove_individuals,
                                   exclude_markers=exclude_markers,
                                   extract_markers=extract_markers,
                                   verbose=verbose)
    }
    if (do.evaluate_check_relatedness) {
      if (verbose) message("Identification of related individuals")
      fail_relatedness <- evaluate_check_relatedness_CRP(qcdir=qcdir,
                                                         name=name,
                                                         imissTh=imissTh,
                                                         highIBDTh=highIBDTh,
                                                         legend_text_size =
                                                           legend_text_size,
                                                         legend_title_size =
                                                           legend_title_size,
                                                         axis_text_size =
                                                           axis_text_size,
                                                         axis_title_size =
                                                           axis_title_size,
                                                         title_size =
                                                           title_size,
                                                         interactive=FALSE)
      write.table(fail_relatedness$failIDs,
                  file=paste(out, ".fail-IBD.IDs", sep=""),
                  row.names=FALSE, quote=FALSE, col.names=FALSE,
                  sep="\t")
      if (!is.null(fail_relatedness$failIDs)  &&
          nrow(fail_relatedness$failIDs) != 0) {
        highIBD <- select(fail_relatedness$failIDs,
                          .data$FID, .data$IID)
      }
      p_relatedness <- fail_relatedness$p_IBD
    }
  }
  if (!dont.check_ancestry) {
    if (do.run_check_ancestry) {
      run <- run_check_ancestry(qcdir=qcdir, indir=indir,
                                prefixMergedDataset=prefixMergedDataset,
                                path2plink=path2plink,
                                showPlinkOutput=showPlinkOutput,
                                keep_individuals=keep_individuals,
                                remove_individuals=remove_individuals,
                                exclude_markers=exclude_markers,
                                extract_markers=extract_markers,
                                verbose=verbose)
    }
    if (do.evaluate_check_ancestry) {
      if (verbose) {
        message("Identification of individuals of divergent ancestry")
      }
      fail_ancestry <- evaluate_check_ancestry(qcdir=qcdir, indir=indir,
                                               name=name,
                                               prefixMergedDataset=
                                                 prefixMergedDataset,
                                               europeanTh=europeanTh,
                                               defaultRefSamples =
                                                 defaultRefSamples,
                                               refSamples=refSamples,
                                               refColors=refColors,
                                               refSamplesFile=
                                                 refSamplesFile,
                                               refColorsFile=
                                                 refColorsFile,
                                               refSamplesIID=
                                                 refSamplesIID,
                                               refSamplesPop=
                                                 refSamplesPop,
                                               refColorsColor=
                                                 refColorsColor,
                                               refColorsPop=
                                                 refColorsPop,
                                               studyColor=studyColor,
                                               highlight_samples =
                                                 highlight_samples,
                                               highlight_type =
                                                 highlight_type,
                                               highlight_text_size =
                                                 highlight_text_size,
                                               highlight_color =
                                                 highlight_color,
                                               highlight_shape =
                                                 highlight_shape,
                                               highlight_legend =
                                                 highlight_legend,
                                               legend_text_size =
                                                 legend_text_size,
                                               legend_title_size =
                                                 legend_title_size,
                                               axis_text_size =
                                                 axis_text_size,
                                               axis_title_size =
                                                 axis_title_size,
                                               title_size = title_size,
                                               interactive=FALSE)
      write.table(fail_ancestry$fail_ancestry,
                  file=paste(out, ".fail-ancestry.IDs",
                             sep=""),
                  quote=FALSE, row.names=FALSE, col.names=FALSE)
      if (!is.null(fail_ancestry$fail_ancestry) &&
          nrow(fail_ancestry$fail_ancestry) != 0) {
        ancestry <- select(fail_ancestry$fail_ancestry,
                           .data$FID, .data$IID)
      }
      p_ancestry <- fail_ancestry$p_ancestry
    }
  }
  
  fail_list <- list(missing_genotype=missing_genotype,
                    highIBD=highIBD,
                    outlying_heterozygosity=outlying_heterozygosity,
                    mismatched_sex=mismatched_sex,
                    ancestry=ancestry)
  
  if(verbose) message(paste("Combine fail IDs into ", out, ".fail.IDs",
                            sep=""))
  
  uniqueFails <- do.call(rbind, fail_list)
  uniqueFails <- uniqueFails[!duplicated(uniqueFails$IID),]
  
  write.table(uniqueFails, file=paste(out, ".fail.IDs",sep=""),
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  plots_sampleQC <- list(p_sexcheck=p_sexcheck,
                         p_het_imiss=p_het_imiss,
                         p_relatedness=p_relatedness,
                         p_ancestry=p_ancestry)
  plots_sampleQC <- plots_sampleQC[sapply(plots_sampleQC,
                                          function(x) !is.null(x))]
  subplotLabels <- LETTERS[1:length(plots_sampleQC)]
  
  if (!is.null(p_ancestry)) {
    ancestry_legend <- cowplot::get_legend(p_ancestry)
    plots_sampleQC$p_ancestry <-  plots_sampleQC$p_ancestry +
      theme(legend.position = "None")
    plots_sampleQC$ancestry_legend <- ancestry_legend
    subplotLabels <- c(subplotLabels, "")
  }
  
  if (!is.null(p_sexcheck) && !is.null(p_het_imiss)) {
    first_plots <- cowplot::plot_grid(plotlist=plots_sampleQC[1:2],
                                      nrow=2,
                                      align = "v",
                                      axis = "lr",
                                      labels=subplotLabels[1:2],
                                      label_size = subplot_label_size
    )
    if (!is.null(p_ancestry)) {
      if (!is.null(p_relatedness)) {
        rel_heights <- c(2, 1, 1, 0.3)
        plots_sampleQC <- list(first_plots,
                               plots_sampleQC[[3]],
                               plots_sampleQC[[4]],
                               plots_sampleQC[[5]]
        )
        subplotLabels <- c("", subplotLabels[3:5])
      } else {
        rel_heights <- c(2, 1, 0.3)
        plots_sampleQC <- list(first_plots,
                               plots_sampleQC[[3]],
                               plots_sampleQC[[4]]
        )
        subplotLabels <- c("", subplotLabels[3:4])
      }
    } else {
      if (!is.null(p_relatedness)) {
        rel_heights <- c(2, 1)
        plots_sampleQC <- list(first_plots, plots_sampleQC[[3]])
        subplotLabels <- c("", subplotLabels[3])
      } else {
        rel_heights <- 1
        plots_sampleQC <- first_plots
        subplotLabels <- ""
      }
    }
  } else {
    if (!is.null(p_ancestry)) {
      rel_heights <- c(rep(1, length(plots_sampleQC) -1), 0.3)
    } else {
      rel_heights <- c(rep(1, length(plots_sampleQC)))
    }
  }
  # p_sampleQC <- cowplot::plot_grid(plotlist=plots_sampleQC,
  #                                  nrow=length(plots_sampleQC),
  #                                  labels=subplotLabels,
  #                                  label_size = subplot_label_size,
  #                                  rel_heights=rel_heights)
  if (interactive) {
    # print(p_sampleQC)
  }
  return(list(fail_list=fail_list #,p_sampleQC=p_sampleQC
              ))
}

checkFormat <- function(prefix) {
  if (!file.exists(paste(prefix, ".fam", sep=""))){
    stop("plink family file: ", prefix, ".fam does not exist.")
  }
  if (!file.exists(paste(prefix, ".bim", sep=""))){
    stop("plink snp file: ", prefix, ".bim does not exist.")
  }
  if (!file.exists(paste(prefix, ".bed", sep=""))){
    stop("plink binary file: ", prefix, ".bed does not exist.")
  }
}

checkRemoveIDs <- function(prefix, remove_individuals=NULL, keep_individuals) {
  removeIDs <- NULL
  allIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                              data.table=FALSE, stringsAsFactors=FALSE,
                              header=FALSE)
  allIDs <- allIDs[,1:2]
  
  if (!is.null(remove_individuals)) {
    if (!file.exists(remove_individuals)) {
      stop("File with individuals to remove from analysis does not exist: ",
           remove_individuals)
    }
    removeIDs <- data.table::fread(remove_individuals, data.table=FALSE,
                                   stringsAsFactors=FALSE,
                                   header=FALSE)
  }
  if (!is.null(keep_individuals)) {
    if (!file.exists(keep_individuals)) {
      stop("File with individuals to keep in analysis does not exist: ",
           keep_individuals)
    }
    pre_keepIDs <- data.table::fread(keep_individuals, data.table=FALSE,
                                     stringsAsFactors=FALSE,
                                     header=FALSE)
    if(ncol(pre_keepIDs) != 2) {
      stop("File keep_individual is not in the right format; should be ",
           "two columns separated by space/tab.")
    }
    removeIDs <- rbind(removeIDs,
                       allIDs[!allIDs[,2] %in% pre_keepIDs[,2],])
  }
  return(removeIDs)
}

cleanData_HS = function (indir, name, qcdir = indir, filterSex = TRUE, filterHeterozygosity = TRUE, 
          filterSampleMissingness = TRUE, filterAncestry = TRUE, filterRelated = TRUE, 
          filterSNPMissingness = TRUE, lmissTh = 0.01, filterHWE = TRUE, 
          hweTh = 1e-05, filterMAF = TRUE, macTh = 20, mafTh = NULL, 
          path2plink = NULL, verbose = FALSE, keep_individuals = NULL, 
          remove_individuals = NULL, exclude_markers = NULL, extract_markers = NULL, 
          showPlinkOutput = TRUE) {
  sampleFilter <- c(filterAncestry, filterRelated, filterSex, 
                    filterHeterozygosity, filterSampleMissingness)
  markerFilter <- c(filterHWE, filterMAF, filterSNPMissingness)
  prefix <- makepath(indir, name)
  out <- makepath(qcdir, name)
  if (!any(c(sampleFilter, markerFilter))) {
    stop("No per-sample and per-marker filters chosen")
  }
  if (!any(sampleFilter)) {
    message("No per-sample filter chosen, carry on with removing markers ", 
            "that fail per-marker QCs")
  }
  if (!any(markerFilter)) {
    message("No per-marker filter chosen, carry on with removing samples ", 
            "that fail per-samples QCs")
  }
  checkFormat(prefix)
  path2plink <- checkPlink(path2plink)
  args_filter <- checkFiltering(extract_markers = extract_markers, 
                                exclude_markers = exclude_markers)
  removeIDs <- checkRemoveIDs(prefix = prefix, remove_individuals = remove_individuals, 
                              keep_individuals = keep_individuals)
  allIDs <- data.table::fread(paste(prefix, ".fam", sep = ""), 
                              data.table = FALSE, stringsAsFactors = FALSE, header = FALSE)
  allIDs <- allIDs[, 1:2]
  if (any(sampleFilter)) {
    if (filterRelated) {
      fail_ibd_ids <- paste0(out, ".fail-IBD.IDs")
      if (!file.exists(fail_ibd_ids)) {
        stop("filterRelated is TRUE but file ", out, 
             ".fail-IBD.IDs does not exist")
      }
      else {
        if (verbose) {
          message("Read individual IDs that failed relatedness check")
        }
        if (file.size(fail_ibd_ids) == 0) {
          if (verbose) {
            message("No individuals failed relatedness check")
          }
        }
        else {
          removeIDs <- rbind(removeIDs, data.table::fread(fail_ibd_ids, 
                                                          data.table = FALSE, stringsAsFactors = FALSE, 
                                                          header = FALSE))
        }
      }
    }
    if (filterAncestry) {
      fail_ancestry_ids <- paste0(out, ".fail-ancestry.IDs")
      if (!file.exists(fail_ancestry_ids)) {
        stop("filterAncestry is TRUE but file ", out, 
             ".fail-ancestry.IDs does not exist")
      }
      else {
        if (verbose) {
          message("Read individual IDs that failed ancestry check")
        }
        if (file.size(fail_ancestry_ids) == 0) {
          if (verbose) {
            message("No individuals failed ancestry check")
          }
        }
        else {
          removeIDs <- rbind(removeIDs, data.table::fread(fail_ancestry_ids, 
                                                          data.table = FALSE, stringsAsFactors = FALSE, 
                                                          header = FALSE))
        }
      }
    }
    if (filterHeterozygosity) {
      fail_het_ids <- paste0(out, ".fail-het.IDs")
      if (!file.exists(fail_het_ids)) {
        stop("filterHeterozygosity is TRUE but file ", 
             out, ".fail-het.IDs does not exist")
      }
      else {
        if (verbose) {
          message("Read individual IDs that failed heterozygosity ", 
                  "check")
        }
        if (file.size(fail_het_ids) == 0) {
          if (verbose) {
            message("No individuals failed heterozygosity check")
          }
        }
        else {
          removeIDs <- rbind(removeIDs, data.table::fread(fail_het_ids, 
                                                          data.table = FALSE, stringsAsFactors = FALSE, 
                                                          header = FALSE))
        }
      }
    }
    if (filterSampleMissingness) {
      fail_imiss_ids <- paste0(out, ".fail-imiss.IDs")
      if (!file.exists(fail_imiss_ids)) {
        stop("filterSampleMissingness is TRUE but file ", 
             out, ".fail-imiss.IDs does not exist")
      }
      else {
        if (verbose) {
          message("Read individual IDs that failed missingness check")
        }
        if (file.size(fail_imiss_ids) == 0) {
          if (verbose) {
            message("No individuals failed missingness check")
          }
        }
        else {
          removeIDs <- rbind(removeIDs, data.table::fread(fail_imiss_ids, 
                                                          data.table = FALSE, stringsAsFactors = FALSE, 
                                                          header = FALSE))
        }
      }
    }
    if (filterSex) {
      fail_sexcheck_ids <- paste0(out, ".fail-sexcheck.IDs")
      if (!file.exists(fail_sexcheck_ids)) {
        stop("filterSex is TRUE but file ", out, ".fail-sexcheck.IDs does not exist")
      }
      else {
        if (verbose) {
          message("Read individual IDs that failed sex check")
        }
        if (file.size(fail_sexcheck_ids) == 0) {
          if (verbose) {
            message("No individuals failed sex check")
          }
        }
        else {
          removeIDs <- rbind(removeIDs, data.table::fread(fail_sexcheck_ids, 
                                                          data.table = FALSE, stringsAsFactors = FALSE, 
                                                          header = FALSE))
        }
      }
    }
    if (!is.null(removeIDs)) {
      removeIDs <- removeIDs[!duplicated(removeIDs), ]
      if (nrow(removeIDs) == nrow(allIDs)) {
        stop("All samples are flagged as .fail.IDs ", 
             "no samples remaining to generate the QCed dataset.")
      }
      if (verbose) 
        message("Write file with remove IDs")
      write.table(removeIDs, paste(out, ".remove.IDs", 
                                   sep = ""), col.names = FALSE, row.names = FALSE, 
                  quote = FALSE)
      remove <- c("--remove", paste(out, ".remove.IDs", 
                                    sep = ""))
      fail_samples <- nrow(removeIDs)
    }
    else {
      remove <- NULL
      fail_samples <- 0
    }
  }
  else {
    remove <- NULL
    fail_samples <- 0
  }
  hwe <- NULL
  maf <- NULL
  missing <- NULL
  if (filterHWE) {
    if (is.null(hweTh)) {
      stop("filterHWE is TRUE but hweTh not specified")
    }
    else {
      hwe <- c("--hwe midp", hweTh)
    }
  }
  if (filterMAF) {
    if (is.null(mafTh) && is.null(macTh)) {
      stop("filterMAF is TRUE but neither mafTh or macTh are provided")
    }
    all_samples <- R.utils::countLines(paste(prefix, ".fam", 
                                             sep = ""))
    keep_samples <- as.numeric(all_samples) - fail_samples
    if (!is.null(mafTh) && !is.null(macTh)) {
      if (verbose) {
        message("Both mafTh and macTh provided, macTh=", 
                macTh, " is used (corresponds to mafTh=", 
                round(mafTh, 6), ")")
      }
    }
    else if (!is.null(mafTh)) {
      if (is.null(macTh)) 
        macTh <- mafTh * (2 * keep_samples)
      if (verbose) {
        message("The mafTh is ", mafTh, " which corresponds to a mcfTh=", 
                macTh)
      }
    }
    else {
      if (is.null(mafTh)) 
        mafTh <- macTh/(2 * keep_samples)
      if (verbose) {
        message("The macTh is ", macTh, " which corresponds to a mafTh=", 
                round(mafTh, 6))
      }
    }
    maf <- c("--maf", mafTh)
  }
  if (filterSNPMissingness) {
    if (is.null(lmissTh)) {
      stop("filterSNPMissingness is TRUE but lmissTh not specified")
    }
    else {
      missing <- c("--geno", lmissTh)
    }
  }
  if (verbose) 
    message("Remove individual IDs and markers IDs that failed QC")
  sys::exec_wait(path2plink, args = c("--bfile", prefix, remove, 
                                      maf, hwe, missing,"--remove", remove_individuals, "--make-bed", "--out", paste(out, 
                                                                                      ".clean", sep = ""), args_filter), std_out = showPlinkOutput, 
                 std_err = showPlinkOutput)
  keepIDs <- data.table::fread(paste0(out, ".clean.fam"), 
                               data.table = FALSE, stringsAsFactors = FALSE, header = FALSE)
  keepIDs <- keepIDs[, 1:2]
  colnames(keepIDs) <- c("FID", "IID")
  colnames(removeIDs) <- c("FID", "IID")
  return(list(passIDs = keepIDs, failIDs = removeIDs))
}

####### 

exclude_relatedness = check_relatedness_CRP(path2plink=path2plink,
                                        indir = indir,
                                        qcdir = qcdir, 
                                        name = name, highIBDTh=0.2,
                           imissTh=0.03, remove_individuals = "fail_IDs.txt")



exclude_related = exclude_relatedness[["plot_data"]] %>% 
  as_tibble() %>% 
  filter(PI_HAT > .2,
         IID1 != IID2) %>% 
  arrange(desc(PI_HAT))





fail_sex <- check_sex(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
                      verbose=TRUE, path2plink=path2plink,
                      maleTh = 0.8, 
                      femaleTh = .2, remove_individuals = "fail_IDs.txt") 
# ggsave(plot = last_plot(), filename = "sex_check.png", width = 6)

fail_sex[["fail_sex"]][["FID"]][fail_sex[["fail_sex"]][["SNPSEX"]]=="0"]

fail_sex[["fail_sex"]] %>% 
  as_tibble() %>% 
  filter(PEDSEX!=SNPSEX,
         SNPSEX != 0)
## warfer029 PEDSEX is wrong. should be 2. 
## warfer075,076,089 have mismatches that won't be filtered to call rate (90%)



sample_ids  <- cleanData(
  indir = indir,
  name = name,
  qcdir = indir,
  filterSex = F,
  filterHeterozygosity = F,
  filterSampleMissingness = F,
  filterAncestry = F,
  filterRelated = F,
  filterSNPMissingness = T,
  lmissTh = .1,
  filterHWE = F,
  filterMAF = F,
  path2plink = path2plink,
  verbose = FALSE,
  keep_individuals = NULL,
  remove_individuals = "fail_IDs.txt",
  exclude_markers = NULL,
  extract_markers = NULL,
  showPlinkOutput = TRUE
)


