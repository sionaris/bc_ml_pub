# The Gene Symbol to Entrez annotation function #####

# It depends on properly defining an ID_Map and Aliases data frame
# on the main script for preprocessing.

annot = function(expressionDF, ID_Map, Aliases){
  # All exprs datasets have "X" as their first column 
  expressionDF = expressionDF %>% dplyr::rename(probe = X)
  
  Join = ID_Map %>% right_join(expressionDF) %>% 
    dplyr::select(probe, EntrezGene.ID, HGNC_Official) %>%
    dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>% distinct()
  Join$EntrezGene.ID = as.numeric(Join$EntrezGene.ID)
  Join = Join[order(Join$EntrezGene.ID),]
  
  dupsg = Join[which(duplicated(Join$EntrezGene.ID)),]
  dupsg_names = as.character(dupsg$EntrezGene.ID)
  dups_entrez = Join[Join$EntrezGene.ID %in% dupsg_names,]
  dups_entrez$EntrezGene.ID = as.numeric(dups_entrez$EntrezGene.ID)
  dups_entrez = dups_entrez[order(dups_entrez$EntrezGene.ID),]
  
  dupsp = Join[which(duplicated(Join$probe)),]
  dupsp_names = as.character(dupsp$probe)
  dups_probes = Join[Join$probe %in% dupsp_names,]
  dups_probes$EntrezGene.ID = as.numeric(dups_probes$EntrezGene.ID)
  dups_probes = dups_probes[order(dups_probes$probe),]
  
  # The full initial set of duplicates is:
  
  full_dups = rbind(dups_entrez, dups_probes)
  full_dups$EntrezGene.ID = as.numeric(full_dups$EntrezGene.ID)
  full_dups = full_dups[order(full_dups$EntrezGene.ID, full_dups$probe),]
  full_dups = full_dups %>% left_join(Aliases, by ="probe") %>% distinct()
  full_dups$final_choice = NA
  
  for(i in 1:nrow(full_dups)){
    if(is.na(full_dups$HGNC_Symbol[i]) == T && full_dups$HGNC_Official[i] == "Yes"){
      full_dups$HGNC_Symbol[i] = full_dups$probe[i]
      full_dups$Entrez[i] = full_dups$EntrezGene.ID[i]
    }
  }
  
  full_dups = full_dups %>% dplyr::filter(is.na(HGNC_Symbol)==F)
  
  # Dealing with duplicates using a complex for loop the philosophy
  # of which is: 
  
  # 1) If official HGNC match and same Entrez ID: 
  # put the probe in final_choice
  
  # 2) If non-official HGNC match but there is a yes somewhere for this 
  # probe: put NA
  
  # 3) If the probe was not an HGNC symbol but an alias
  # check if it has been considered anywhere as an HGNC symbol, if not
  # then check if there's only one probe matched to this HGNC symbol
  # and if yes check that this is the only probe name (redundant). If 
  # all are true final_choice is that HGNC symbol
  
  # 4) If the probe was not an HGNC symbol but an alias
  # check if it has been considered anywhere as an HGNC symbol, if not
  # then check if there's only one HGNC symbol matched to this probe
  # and if yes check if the HGNC symbol's name is that name (redundant).
  # If all are true then put that HGNC symbol as final_choice
  
  # ELSE put NA
  for (i in 1:nrow(full_dups)){
    if (full_dups$HGNC_Official[i] == "Yes" &&
        full_dups$EntrezGene.ID[i] == full_dups$Entrez[i]){
      full_dups$final_choice[i] = full_dups$probe[i]
    } else if (full_dups$HGNC_Official[i] == "No" &&
               length(unique(full_dups$HGNC_Official[full_dups$probe == full_dups$probe[i]])) > 1){
      full_dups$final_choice[i] = NA
    } else if (full_dups$HGNC_Official[i] == "No" &&
               length(unique(full_dups$HGNC_Official[full_dups$probe == full_dups$probe[i]])) == 1 &&
               length(unique(full_dups$probe[full_dups$HGNC_Symbol == full_dups$HGNC_Symbol[i]])) == 1 &&
               full_dups$EntrezGene.ID[i] == full_dups$Entrez[i]){
      full_dups$final_choice[i] = full_dups$HGNC_Symbol[i]
    } else if (full_dups$HGNC_Official[i] == "No" &&
               length(unique(full_dups$HGNC_Official[full_dups$probe == full_dups$probe[i]])) == 1 &&
               length(unique(full_dups$HGNC_Symbol[full_dups$probe == full_dups$probe[i]])) == 1 &&
               length(unique(full_dups$probe[full_dups$EntrezGene.ID == full_dups$EntrezGene.ID[i]])) == 1){
      full_dups$final_choice[i] = full_dups$HGNC_Symbol[i]
    } else {
      full_dups$final_choice[i] = NA
    }
  }
  
  # Check for unmatched probes
  nonas_names = as.character(unique(full_dups$probe[which(is.na(full_dups$final_choice) == FALSE)]))
  nas_names = as.character(unique(full_dups$probe[which(is.na(full_dups$final_choice) == TRUE)]))
  check = nas_names %in% nonas_names
  Unmatched1 = as.data.frame(nas_names[which(check == FALSE)])
  colnames(Unmatched1) = "Unmatched probes"
  
  # Ignore the final_choice values of probes that have been matched
  # elsewhere
  for (b in 1:nrow(full_dups)){
    if (full_dups$probe[b] %in% nonas_names &&
        is.na(full_dups$final_choice[b]) == TRUE){
      full_dups$final_choice[b] = "Ignore"
    }
  }
  
  # Same but using HGNC symbols so that we ignore the HGNC symbols
  # that have been matched
  for (c in 1:length(full_dups$probe)){
    if (full_dups$HGNC_Symbol[c] %in% nonas_names &&
        is.na(full_dups$final_choice[c]) == TRUE){
      full_dups$final_choice[c] = "Ignore"
    }
  }
  
  # Find all probes which do not have final_choice and match them to 
  # the corresponding HGNC symbol if this is unique and not present
  # in final choice. Create a match only when EntrezGene.ID and Entrez
  # are the same
  for (e in 1:length(full_dups$probe)){
    if (is.na(full_dups$final_choice[e]) == TRUE &&
        full_dups$HGNC_Official[e] == "No" &&
        length(unique(full_dups$HGNC_Official[full_dups$probe == full_dups$probe[e]])) == 1 &&
        full_dups$EntrezGene.ID[e] == full_dups$Entrez[e] &&
        nrow(distinct(full_dups[full_dups$probe==full_dups$probe[e] &
                                is.na(full_dups$final_choice) == TRUE,c(1,4)])) == 1){
      full_dups$final_choice[e] = full_dups$HGNC_Symbol[e]
    }
  }
  
  # Check for unmatched probes vol.2
  nonas_names = as.character(unique(full_dups$probe[which(is.na(full_dups$final_choice) == FALSE)]))
  nas_names = as.character(unique(full_dups$probe[which(is.na(full_dups$final_choice) == TRUE)]))
  check = nas_names %in% nonas_names
  Unmatched1 = as.data.frame(nas_names[which(check == FALSE)])
  colnames(Unmatched1) = "Unmatched probes"
  full_dups$final_choice[full_dups$probe %in% Unmatched1$`Unmatched probes`] = "Unclear"
  
  # The most challenging of probes
  rogues = full_dups %>% dplyr::filter(final_choice == "Ignore" || final_choice == "Unclear")
  rogue_names = unique(rogues$probe)
  rights = full_dups %>% dplyr::filter(is.na(final_choice) == FALSE) %>%
    dplyr::filter(final_choice != "Ignore") %>%
    dplyr::filter(final_choice != "Unclear")
  right_names = unique(rights$probe)
  real_rogues = rogue_names[which(rogue_names %in% right_names == FALSE)]
  bandits = unique(full_dups$probe[full_dups$probe %in% real_rogues]) %>%
    as.data.frame()
  final_check = Join[Join$probe %in% bandits$.,]
  full_dups$EntrezGene.ID[full_dups$probe %in% final_check$probe] = NA
  
  # Filter the full_dups dataset in order to use it to 
  # correct the Join dataset
  full_dups_filt = full_dups %>% dplyr::filter(is.na(final_choice) == FALSE) %>%
    dplyr::filter(final_choice != "Unclear") %>%
    dplyr::filter(final_choice != "Ignore") %>%
    dplyr::select(probe, Entrez, HGNC_Official, final_choice) %>% distinct()
  
  Unmatched2 = as.data.frame(full_dups_filt$probe[which(duplicated(full_dups_filt$probe))])
  colnames(Unmatched2) = "Unmatched probes"
  Unmatched = unique(rbind(Unmatched1, Unmatched2))
  Unmatched_names = unique(Unmatched$`Unmatched probes`)
  
  dup_positions = which(duplicated(full_dups_filt$probe))
  dup_names = unique(full_dups_filt$probe[dup_positions])
  full_dups_filt$Entrez[full_dups_filt$probe %in% dup_names] = NA
  full_dups_filt$Entrez[full_dups_filt$probe %in% final_check$probe] = NA
  Join$EntrezGene.ID[Join$probe %in% unique(full_dups$probe)] = "Remove"
  Join$EntrezGene.ID[Join$probe %in% Unmatched_names] = "Remove"
  Join_pre = Join %>% dplyr::filter(EntrezGene.ID != "Remove") %>%
    dplyr::select(-HGNC_Official)
  
  # HGNC to HGNC annotation preparations for intClust
  hgnc_annot_filt = full_dups_filt %>% dplyr::select(probe, final_choice) %>%
    dplyr::rename(Gene.Symbol = final_choice)
  hgnc_annot_pre = Join_pre %>% dplyr::select(probe)
  hgnc_annot_pre$Gene.Symbol = as.character(hgnc_annot_pre$probe)
  hgnc_annot_pre2 = rbind(hgnc_annot_pre, hgnc_annot_filt)
  uniq2 = unique(hgnc_annot_pre2)
  hgnc_annot_pre2$Gene.Symbol[hgnc_annot_pre2$probe %in% uniq2$probe[which(duplicated(uniq2$probe) == TRUE)]] = NA
  hgnc_annot_pre2 = unique(hgnc_annot_pre2)
  hgnc_annot = as.data.frame(expressionDF[,1])
  colnames(hgnc_annot) = "probe"
  hgnc_annot = hgnc_annot %>% left_join(hgnc_annot_pre2)
  
  # The following is added because CDCA1 and NUF2 create confusions
  # in most datasets especially in PAM50
  if ("CDCA1" %in% hgnc_annot$probe &&
      "NUF2" %in% hgnc_annot$probe){
    hgnc_annot$Gene.Symbol[hgnc_annot$probe == "CDCA1"] = NA
  }
  
  # Final steps for annotation data frames
  full_dups_filt = full_dups_filt %>% dplyr::select(probe, Entrez) %>%
    dplyr::rename(EntrezGene.ID=Entrez)
  Join_final = unique(rbind(Join_pre, full_dups_filt))
  Join_annot = as.data.frame(expressionDF[,1])
  colnames(Join_annot) = "probe"
  Join_annot2 = Join_annot %>% left_join(Join_final)
  Join_annot2$EntrezGene.ID[Join_annot2$probe %in% final_check$probe] = NA
  Join_annot2$EntrezGene.ID[Join_annot2$probe %in% dup_names] = NA
  dannot = distinct(Join_annot2)
  
  # Do a PAM50 check
  data("pam50")
  inconsistencies = as.data.frame(pam50$centroids.map[,c("probe", "EntrezGene.ID")])
  inconsistencies[51,] = c("AURKA", "6790")
  rownames(inconsistencies) = as.character(inconsistencies$probe)
  dannot$f_probe = as.character(dannot$probe)
  for (j in 1:length(inconsistencies$EntrezGene.ID)){
    dannot$f_probe[dannot$EntrezGene.ID == inconsistencies$EntrezGene.ID[j]] = inconsistencies$probe[j]
  }
  
  if ("CDCA1" %in% dannot$probe && 
      "NUF2" %in% dannot$probe){
    dannot$EntrezGene.ID[dannot$probe == "CDCA1"] = NA
  }
  
  sigs = length(intersect(dannot$EntrezGene.ID, pam50$centroids.map$EntrezGene.ID))
  message(paste0(sigs, " pam50 signatures were found"))
  
  dannot$EntrezGene.ID[dannot$probe %in% final_check$probe] = NA
  dannot$EntrezGene.ID[dannot$probe %in% dup_names] = NA
  
  expressionDF = expressionDF %>% inner_join(dannot, by = "probe") 
  
  # If after all this process there still are duplicate probes
  # keep the most variant probe for this gene
  last_dups = dannot$f_probe[which(duplicated(dannot$f_probe) == TRUE)]
  last_dups_probe = dannot$probe[dannot$f_probe %in% unique(last_dups)]
  if (length(last_dups) != 0){
    means = expressionDF[expressionDF$probe %in% last_dups_probe,]
    means = means %>% dplyr::select(-probe, -EntrezGene.ID) %>% group_by(f_probe) %>%
      summarise_all(mean, na.rm = TRUE)
    for(i in 1:nrow(expressionDF)){
      if(expressionDF$f_probe[i] %in% means$f_probe){
        expressionDF[i,2:(ncol(expressionDF)-2)] = means[means$f_probe==expressionDF$f_probe[i],
                                                         2:ncol(means)]
      }
    }
    expressionDF = expressionDF %>%    
      dplyr::select(-f_probe, -EntrezGene.ID) %>% inner_join(dannot)
    if ("CDCA1" %in% dannot$probe && 
        "NUF2" %in% dannot$probe){
      expressionDF$f_probe[expressionDF$probe == "CDCA1"] = "CDCA_ALIAS"
    }
    expressionDF = expressionDF %>%
      dplyr::select(-probe, -EntrezGene.ID) %>%
      dplyr::rename(probe=f_probe) %>%
      dplyr::select(probe, everything())
  } else {
    expressionDF = expressionDF %>%
      dplyr::select(-probe, -EntrezGene.ID) %>%
      dplyr::rename(probe=f_probe) %>%
      dplyr::select(probe, everything())
  }
  
  if ("CDCA1" %in% dannot$probe && 
      "NUF2" %in% dannot$probe){
    dannot$f_probe[dannot$probe == "CDCA1"] = "CDCA_ALIAS"
  }
  dannot = dannot %>% dplyr::select(-probe) %>%
    dplyr::rename(probe=f_probe) %>% distinct()
  rownames(dannot) = as.character(dannot[,"probe"])
  
  message("Annotation data frame created successfully")
  Entrez = expressionDF %>% inner_join(dannot) %>%
    dplyr::select(-probe) %>%
    dplyr::select(EntrezGene.ID, everything())
  HGNC = expressionDF %>% inner_join(hgnc_annot) %>%
    dplyr::select(-probe) %>%
    dplyr::select(Gene.Symbol, everything())
  return(list(Entrez = as.data.frame(Entrez), HGNC = as.data.frame(HGNC), 
              dannot_Entrez = as.data.frame(dannot), 
              dannot_HGNC = as.data.frame(hgnc_annot)))
}

# The Genefu predictions function - full #####
Genefu_predictions <- function(expressionDF, HGNC, dannot_Entrez, dannot_HGNC){
  library(genefu)
  library(dplyr)
  library(org.Hs.eg.db)
  
  # The purpose of this function is to generate genefu predictions
  # based on properly Entrez-annotated datasets. expressionDF is the Entrez
  # exprs and HGNC is the HGNC exprs for intClust
  
  entrezgenes = as.character(expressionDF$EntrezGene.ID)
  hgncgenes = as.character(HGNC$Gene.Symbol)
  
  expressionDF2 <- expressionDF %>% dplyr::select(-EntrezGene.ID)
  HGNC2 <- HGNC %>% dplyr::select(-Gene.Symbol)
  
  exprsDF <- t(expressionDF2)
  exprsHGNC <- t(HGNC2)
  
  colnames(exprsDF) = entrezgenes
  colnames(exprsHGNC) = hgncgenes
  
  # Predictions
  # Mammaprint
  # std = "robust": Standardization based on the 0.025 and 0.975 quantiles
  # same as in the case of 
  data("sig.gene70")
  RNGversion("4.0.2")
  set.seed(123)
  mammaprint <- gene70(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE,
                       std = "robust")
  mammaprint_score <- mammaprint$score %>% as.data.frame()
  mammaprint_risk <- mammaprint$risk %>% as.data.frame()
  message("Mammaprint: OK")
  
  # rorS
  RNGversion("4.0.2")
  set.seed(123)
  RORS <- rorS(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE)
  rorS_score <- RORS$score %>% as.data.frame()
  rorS_risk <- RORS$risk %>% as.data.frame()
  message("rorS: OK")
  
  # EndoPredict
  # This algorithm performs the following rescaling process for values to 
  # approximate the scale of the original Affymetrix data:
  # https://rdrr.io/bioc/genefu/src/R/endoPredict.R
  
  data("sig.endoPredict")
  RNGversion("4.0.2")
  set.seed(123)
  EndoPredict = endoPredict(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE)
  endoPredict_score <- EndoPredict$score %>% as.data.frame()
  endoPredict_risk <- EndoPredict$risk %>% as.data.frame()
  message("EndoPredict: OK")
  
  # Oncotype DX
  data("sig.oncotypedx")
  RNGversion("4.0.2")
  set.seed(123)
  OncotypeDX = oncotypedx(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE,
                          do.scaling = TRUE)
  oncotypedx_score <- OncotypeDX$score %>% as.data.frame()
  oncotypedx_risk <- OncotypeDX$risk %>% as.data.frame()
  message("Oncotype DX: OK")
  
  # When choosing ".robust" versions of the algorithm, scaling of the expression
  # values is performed using the rescale() function of the genefu package. This function
  # rescales values x based on quantiles specified by the user such that 
  # x’ = (x - q1) / (q2 - q1) where q is the specified quantile, q1 = q / 2, 
  # q2 = 1 - q/2) and x’ are the new rescaled values. This version of scaling is
  # more concordant with clinical experience. The robust version of the algorithm
  # is by default chosen when these algorithms are run.
  
  # SCMGENE
  data("scmgene.robust")
  scmgene <- molecular.subtyping(sbt.model = "scmgene", data = exprsDF, 
                                 annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("scmgene: OK")
  
  # SCMOD1
  data("scmod1.robust")
  scmod1 <- molecular.subtyping(sbt.model = "scmod1", data = exprsDF,
                                annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("scmod1: OK")
  
  # SCMOD2
  data("scmod2.robust")
  scmod2 <- molecular.subtyping(sbt.model = "scmod2", data = exprsDF,
                                annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("scmod2: OK")
  
  # PAM50
  data("pam50.robust")
  PAM50 <- molecular.subtyping(sbt.model = "pam50", data = exprsDF,
                               annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("pam50: OK")
  
  # ClaudinLow
  # This algorithm by default performs median centering
  data("claudinLowData")
  ClaudinLow <- molecular.subtyping(sbt.model = "claudinLow", data = exprsDF,
                                    annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("ClaudinLow: OK")
  
  # SSP2006
  data("ssp2006.robust")
  ssp_2006 <- molecular.subtyping(sbt.model = "ssp2006", data = exprsDF,
                                  annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("ssp2006: OK")
  
  # SSP2003
  data("ssp2003.robust")
  ssp_2003 <- molecular.subtyping(sbt.model = "ssp2003", data = exprsDF,
                                  annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("ssp2003: OK")
  
  # GENIUS
  # rescaling here is performed to yield more robust comparisons with AURKA
  # ESR1 and ERBB2 modules
  data("sig.genius")
  RNGversion("4.0.2")
  set.seed(123)
  GENIUS_score <- genius(data = exprsDF, annot = dannot_Entrez,
                         do.mapping = TRUE, do.scale = TRUE)$score %>% 
    as.data.frame()
  
  message("GENIUS: OK")
  
  # GGI - similar logic to the Nottingham Index score
  data("sig.ggi")
  RNGversion("4.0.2")
  set.seed(123)
  GGI <- ggi(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE)
  GGI_score <- GGI$score %>% as.data.frame()
  GGI_risk <- GGI$risk %>% as.data.frame()
  
  message("GGI: OK")
  
  # PIK3CA_GS
  data("sig.pik3cags")
  RNGversion("4.0.2")
  set.seed(123)
  PIK3CA_GS <- pik3cags(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE) %>% 
    as.data.frame()
  
  message("PIK3CA_GS: OK")
  
  # IC10
  # This algorithm normalises features by default (z-score)
  intClust <- molecular.subtyping(sbt.model = "intClust", data = exprsHGNC,
                                  annot = dannot_HGNC, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("iC10 OK")
  
  # Finalising
  predictions <- cbind(scmgene, scmod1, scmod2, PAM50, ClaudinLow, ssp_2006, 
                       ssp_2003, GENIUS_score, GGI_score, GGI_risk, PIK3CA_GS, 
                       mammaprint_risk, mammaprint_score, rorS_risk, rorS_score, 
                       endoPredict_risk, endoPredict_score, oncotypedx_risk, 
                       oncotypedx_score, intClust)
  colnames(predictions) <- c("scmgene", "scmod1", "scmod2", "pam50", "ClaudinLow",
                             "ssp2006",  "ssp2003", "GENIUS_score",  "GGI_score",
                             "GGI_risk", "PIK3CA_GS", "Mammaprint_risk", 
                             "Mammaprint_score", "rorS_risk", "rorS_score",
                             "EndoPredict_risk", "EndoPredict_score", 
                             "OncotypeDX_risk", "OncotypeDX_score", "IC10")
  predictions$Sample.ID = as.character(rownames(predictions))
  predictions <- predictions %>% dplyr::select(Sample.ID, everything())
  return(predictions)
}
# Feed one expression matrix.

# The Genefu predictions function - fewer predictors #####
Genefu_predictions_2 <- function(expressionDF, HGNC, dannot_Entrez, dannot_HGNC){
  library(genefu)
  library(dplyr)
  library(org.Hs.eg.db)
  
  # The purpose of this function is to generate genefu predictions
  # based on properly Entrez-annotated datasets. expressionDF is the Entrez
  # exprs and HGNC is the HGNC exprs for intClust
  
  entrezgenes = as.character(expressionDF$EntrezGene.ID)
  hgncgenes = as.character(HGNC$Gene.Symbol)
  
  expressionDF2 <- expressionDF %>% dplyr::select(-EntrezGene.ID)
  HGNC2 <- HGNC %>% dplyr::select(-Gene.Symbol)
  
  exprsDF <- t(expressionDF2)
  exprsHGNC <- t(HGNC2)
  
  colnames(exprsDF) = entrezgenes
  colnames(exprsHGNC) = hgncgenes
  
  # Predictions
  # Mammaprint
  # std = "robust": Standardization based on the 0.025 and 0.975 quantiles
  # same as in the case of 
  data("sig.gene70")
  RNGversion("4.0.2")
  set.seed(123)
  mammaprint <- gene70(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE,
                       std = "robust")
  mammaprint_score <- mammaprint$score %>% as.data.frame()
  mammaprint_risk <- mammaprint$risk %>% as.data.frame()
  message("Mammaprint: OK")
  
  # rorS
  RNGversion("4.0.2")
  set.seed(123)
  RORS <- rorS(data = exprsDF, annot = dannot_Entrez, do.mapping = TRUE)
  rorS_score <- RORS$score %>% as.data.frame()
  rorS_risk <- RORS$risk %>% as.data.frame()
  message("rorS: OK")
  
  # When choosing ".robust" versions of the algorithm, scaling of the expression
  # values is performed using the rescale() function of the genefu package. This function
  # rescales values x based on quantiles specified by the user such that 
  # x’ = (x - q1) / (q2 - q1) where q is the specified quantile, q1 = q / 2, 
  # q2 = 1 - q/2) and x’ are the new rescaled values. This version of scaling is
  # more concordant with clinical experience. The robust version of the algorithm
  # is by default chosen when these algorithms are run.
  
  # SCMOD1
  data("scmod1.robust")
  scmod1 <- molecular.subtyping(sbt.model = "scmod1", data = exprsDF,
                                annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("scmod1: OK")
  
  # PAM50
  data("pam50.robust")
  PAM50 <- molecular.subtyping(sbt.model = "pam50", data = exprsDF,
                               annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("pam50: OK")
  
  # ClaudinLow
  # This algorithm by default performs median centering
  data("claudinLowData")
  ClaudinLow <- molecular.subtyping(sbt.model = "claudinLow", data = exprsDF,
                                    annot = dannot_Entrez, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("ClaudinLow: OK")
  
  # IC10
  # This algorithm normalises features by default (z-score)
  intClust <- molecular.subtyping(sbt.model = "intClust", data = exprsHGNC,
                                  annot = dannot_HGNC, do.mapping = TRUE)$subtype %>% 
    as.data.frame()
  
  message("iC10 OK")
  
  # Finalising
  predictions <- cbind(scmod1, PAM50, ClaudinLow, mammaprint_risk, mammaprint_score,
                       rorS_risk, rorS_score, intClust)
  colnames(predictions) <- c("scmod1", "pam50", "ClaudinLow", "Mammaprint_risk", 
                             "Mammaprint_score", "rorS_risk", "rorS_score", "IC10")
  predictions$Sample.ID = as.character(rownames(predictions))
  predictions <- predictions %>% dplyr::select(Sample.ID, everything())
  return(predictions)
}
# Feed one expression matrix.