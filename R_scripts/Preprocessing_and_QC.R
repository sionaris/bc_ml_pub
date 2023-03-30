# Libraries #####
library(limma)
library(org.Hs.eg.db)
library(dplyr)
library(openxlsx)
library(GEOquery)
library(genefu)
library(stringr)
library(tidyr)
library(EnhancedVolcano)
library(readr)
# library(knitr)
# library(rJava)
readr::local_edition(1) # required for v>=2.60.0 of GEOquery

# Loading data #####

# Downloading the GEO objects for all studies of interest
tvt_series = c("GSE32603", "GSE20181", "GSE55374", "GSE59515", "GSE87411",
               "GSE105777", "GSE126870")

ext_val_series = c("GSE18728", "GSE119262", "GSE111563")

GEOsets_pre = list()
for (i in 1:length(tvt_series)){
  GEOsets_pre[[i]] = unlist(getGEO(tvt_series[i])[[1]])
}
names(GEOsets_pre) = tvt_series; rm(tvt_series, i); gc()

Ext_val = list()
for (i in 1:length(ext_val_series)){
  Ext_val[[i]] = unlist(getGEO(ext_val_series[i])[[1]])
}
names(Ext_val) = ext_val_series; rm(ext_val_series, i); gc()

# Expression matrices for the Bownes, Dunbier and Park studies will be loaded 
# as local files. We stored the downloaded objects in "data/Expression Matrices"
# Loading the remaining Expression Matrices
detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)
Bownes = read.csv("data/Expression Matrices/matrixC1.csv") # available in the andysims lab
Dunbier = read.csv("data/Expression Matrices/matrix_E2.csv") # not available online - provided by Anita Dunbier to andysims lab
Park = read.csv("data/Expression Matrices/matrixC3.csv") # GEOquery fails in retrieving the data
# Removing some unidentified probes
Park = Park %>% dplyr::slice(1:16692,16706:30373)

GEOsets = list(Bownes, GEOsets_pre[[1]], Park, GEOsets_pre[[2]], GEOsets_pre[[3]], 
               GEOsets_pre[[4]], Dunbier, GEOsets_pre[[5]], GEOsets_pre[[6]],
               GEOsets_pre[[7]])
names(GEOsets) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                   "E4_1", "E4_2")

rm(GEOsets_pre, Bownes, Park, Dunbier)

# Phenotypic annotation #####
# Entrez_pheno (HGNC_pheno is the same essentially)
Entrez_pheno = GEOsets
for (i in 1:length(GEOsets)){
  Entrez_pheno[[i]] = as.data.frame(colnames(GEOsets[[i]][,2:ncol(GEOsets[[i]])]))
  colnames(Entrez_pheno[[i]]) = "Sample.ID"
  Entrez_pheno[[i]]$Dataset = names(Entrez_pheno[i])
  if(i <= 3){Entrez_pheno[[i]]$Treatment = "Chemotherapy"} else 
  {Entrez_pheno[[i]]$Treatment = "Endocrine_treatment"}
}

# Starting from C1 - info taken from publication and exprs sample names
Entrez_pheno[[1]]$Timepoint_coded = str_sub(Entrez_pheno[[1]]$Sample.ID, -2, -1)
Entrez_pheno[[1]]$Timepoint = as.character(Entrez_pheno[[1]]$Timepoint_coded)
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T2"] = "2 weeks"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T3"] = "1.5 months"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T4"] = "Surgery"
Entrez_pheno[[1]]$Timepoint_coded[Entrez_pheno[[1]]$Timepoint == "1.5 months"] = "T2.5"
Entrez_pheno[[1]]$Timepoint_coded[Entrez_pheno[[1]]$Timepoint == "Surgery"] = "T3"
Entrez_pheno[[1]]$Patient.ID = str_sub(Entrez_pheno[[1]]$Sample.ID, 1,
                                       nchar(Entrez_pheno[[1]]$Sample.ID) - 3)
str_sub(Entrez_pheno[[1]]$Patient.ID, 1, 3) <- "C1."
Entrez_pheno[[1]]$Response_type = "pCR"
Entrez_pheno[[1]]$Response = NA

# C2
# times
pdata2 = pData(GEOsets[["C2"]])
times = data.frame(rownames(pdata2),pdata2$characteristics_ch1.3)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times$Timepoint_coded = str_sub(times$Timepoint_coded, -2, -1)
times[which(times$Timepoint_coded == "TS"),2] = "T3"
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% inner_join(times, by = "Sample.ID")
Entrez_pheno[[2]]$Timepoint = as.character(Entrez_pheno[[2]]$Timepoint_coded)
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T2"] = "24-96 hours"
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T3"] = "Surgery"

# patient ID's
patient_ids = data.frame(rownames(pdata2), pdata2$characteristics_ch1)
patient_ids = patient_ids %>%
  dplyr::rename(Sample.ID=colnames(patient_ids)[1], 
                Patient.ID=colnames(patient_ids)[2])
patient_ids$Patient.ID = str_sub(patient_ids$Patient.ID, -4, -1)
str_sub(patient_ids$Patient.ID, 0, 0) = "C2."
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% left_join(patient_ids, by = "Sample.ID")

# Response
Entrez_pheno[[2]]$Response_type = "RCB, cutoff: >1"
response = data.frame(rownames(pdata2), pdata2$characteristics_ch1.9)
response = response %>%
  dplyr::rename(Sample.ID=colnames(response)[1], 
                Response=colnames(response)[2])
response$Response = str_sub(response$Response, 11, -1)
response$Response = as.numeric(response$Response)
response[which(response$Response == 3),2] = "Non_responder"
response[which(response$Response == 2),2] = "Non_responder"
response[which(response$Response == 1),2] = "Responder"
response[which(response$Response == 0),2] = "Responder"
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% inner_join(response, by = "Sample.ID")

rm(pdata2, times, patient_ids, response)

# C3
# times
library(GEOquery)
readr::local_edition(1)
pdata3 = pData(getGEO('GSE123845',GSEMatrix=TRUE)[[1]])
times = data.frame(pdata3$title, pdata3$characteristics_ch1.54)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times$Timepoint_coded = str_sub(times$Timepoint_coded, -2, -1)
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(times, by = "Sample.ID")
Entrez_pheno[[3]]$Timepoint = as.character(Entrez_pheno[[3]]$Timepoint_coded)
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T2"] = "3 weeks"
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T3"] = "6 months/Surgery"

# patient ID's
patient_ids = data.frame(pdata3$title)
patient_ids$Patient.ID = as.character(patient_ids$pdata3.title)
patient_ids = patient_ids %>%
  dplyr::rename(Sample.ID=colnames(patient_ids)[1])
patient_ids$Patient.ID = str_sub(patient_ids$Patient.ID, 4, -3)
str_sub(patient_ids$Patient.ID, 0, 0) = "C3."
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(patient_ids, by = "Sample.ID")

# Response
response = data.frame(pdata3$title, pdata3$characteristics_ch1.35)
response = response %>%
  dplyr::rename(Sample.ID=colnames(response)[1], 
                Response=colnames(response)[2])
response$Response_type = "pCR/RD"
response[which(response$Response == "pcr_status: 0"),2] = "Responder"
response[which(response$Response == "pcr_status: 1"),2] = "Non_responder"
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(response, by = "Sample.ID")

rm(pdata3, times, patient_ids, response)
detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)

# E1_1
# times & patient ID's & response
pdata4 = pData(GEOsets[["E1_1"]])
times = data.frame(rownames(pdata4), pdata4$title)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times = separate(times, Timepoint_coded, into = c("Timepoint_coded",
                                                  "col2", "col3", "col4",
                                                  "Response"),
                 sep = ";")
times$Response[117:176] = times$col2[117:176]
times = times %>% dplyr::select(-col2, -col3, -col4)
times = separate(times, Timepoint_coded, into = c("Timepoint_coded",
                                                  "col2", "col3", "col4"),
                 sep = ":")
times$Response_type = "US_tumour_volume"
times[which(times$Response == " responder"), 6] = "Responder"
times[which(times$Response == " nonresponder"), 6] = "Non_responder"
times[which(times$Response == " not assessable"), 6] = NA

times$Patient.ID = str_sub(times$Timepoint_coded,0,-2)
str_sub(times$Patient.ID, 0, 0) = "E1.1."
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "A"),2] = "T1"
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "B"),2] = "T2"
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "C"),2] = "T3"
times$Timepoint = as.character(times$Timepoint_coded)
times$Timepoint[times$Timepoint_coded == "T1"] = "Pre-treatment"
times$Timepoint[times$Timepoint_coded == "T2"] = "2 weeks"
times$Timepoint[times$Timepoint_coded == "T3"] = "Surgery"
times = times %>% dplyr::select(Sample.ID, Timepoint_coded, Timepoint,
                                Patient.ID, Response, Response_type)
Entrez_pheno[[4]] = Entrez_pheno[[4]] %>% inner_join(times, by = "Sample.ID")

rm(pdata4, times)

# E1_2
pdata5 = pData(GEOsets[["E1_2"]])
add = pdata5 %>% dplyr::select(geo_accession, `clinical response:ch1`,
                               `subject:ch1`, `treatment:ch1`) %>%
  dplyr::rename(Sample.ID=geo_accession, Response=`clinical response:ch1`,
                Patient.ID=`subject:ch1`, Timepoint=`treatment:ch1`)
add$Response_type = "US_tumour_volume"
str_sub(add$Patient.ID, 0, -4) = "E1.2."
add[which(add$Timepoint == "pre-treatment"), 4] = "Pre-treatment"
add[which(add$Timepoint == "2 wk endocrine therapy"), 4] = "2 weeks"
add[which(add$Timepoint == "3 mo endocrine therapy"), 4] = "3 months/Surgery"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[which(add$Timepoint == "3 months/Surgery")] = "T3"
add$Timepoint_coded[which(add$Timepoint == "2 weeks")] = "T2"
add$Timepoint_coded[which(add$Timepoint == "Pre-treatment")] = "T1"
Entrez_pheno[[5]] = Entrez_pheno[[5]] %>% inner_join(add, by = "Sample.ID")

rm(add, pdata5)

# E1_3
pdata6 = pData(GEOsets[["E1_3"]])
add = pdata6 %>% dplyr::select(geo_accession, `clinical response:ch1`,
                               `patient id:ch1`, `time point:ch1`) %>%
  dplyr::rename(Sample.ID=geo_accession, Response=`clinical response:ch1`,
                Patient.ID=`patient id:ch1`, Timepoint=`time point:ch1`)
add$Response_type = "US_tumour_volume"
str_sub(add$Patient.ID, 0, 0) = "E1.3."
add[which(add$Timepoint == "pre-treatment"), 4] = "Pre-treatment"
add[which(add$Timepoint == "2 wks of letrozole treatment"), 4] = "2 weeks"
add[which(add$Timepoint == "3 months of letrozole treatment"), 4] = "3 months/Surgery"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[which(add$Timepoint == "3 months/Surgery")] = "T3"
add$Timepoint_coded[which(add$Timepoint == "2 weeks")] = "T2"
add$Timepoint_coded[which(add$Timepoint == "Pre-treatment")] = "T1"
Entrez_pheno[[6]] = Entrez_pheno[[6]] %>% inner_join(add, by = "Sample.ID")
Entrez_pheno[[6]]$Response = gsub("Non-responder", "Non_responder", Entrez_pheno[[6]]$Response)

rm(add, pdata6)

# E2
# E2_supp is created in MS Excel from data again available at the Sims lab space
# provided by Anita Dunbier from the Royal Marsden study ("E2")
E2_supp = openxlsx::read.xlsx("data/Supplementary/E2_extras.xlsx", sheet = 1, colNames = T)
E2_supp = E2_supp %>%
  dplyr::select(Sample.ID, Patient.ID, Timepoint,
                `FactorValue[ClinicalInformation-Response]`) %>%
  dplyr::rename(Response = `FactorValue[ClinicalInformation-Response]`)
str_sub(E2_supp$Patient.ID, 0, 0) = "E2."
E2_supp$Timepoint_coded = as.character(E2_supp$Timepoint)
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "16week post-treatment")] = "T3"
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "2 week post-treatment")] = "T2"
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "pre-treatment")] = "T1"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T1")] = "Pre-treatment"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T2")] = "2 weeks"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T3")] = "16 weeks"
E2_supp$Response = gsub("responder", "Responder", E2_supp$Response)
E2_supp$Response = gsub("non-Responder", "Non_responder", E2_supp$Response)
E2_supp$Response_type = "Ki67 expression"
Entrez_pheno[[7]] = Entrez_pheno[[7]] %>% left_join(E2_supp, by = "Sample.ID")

rm(E2_supp)

# E3
# This set contains reanalysed samples from studies not included in our list
# of studies, so reanalysed samples are kept in this case

pdata8 = pData(GEOsets[["E3"]])
add = pdata8 %>% 
  dplyr::select(`geo_accession`, `source_name_ch2`,`studynum:ch2`, 
                `endocrine therapy response group:ch2`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Timepoint = `source_name_ch2`,
                Patient.ID = `studynum:ch2`, 
                Response = `endocrine therapy response group:ch2`)
str_sub(add$Patient.ID, 0, 0) = "E3."
add[which(add$Response == "resistant"), 4] = "Non_responder"
add[which(add$Response == "sensitive"), 4] = "Responder"
add$Response_type = "pCR (+Ki67)"
add = separate(add, Timepoint, into = c("col1", "col2", "col3",
                                        "Timepoint", "col5"),
               sep = "_")
add = add %>% dplyr::select(-col1, -col2, -col3, -col5)
add[which(add$Timepoint == "BL"), 2] = "Pre-treatment"
add[which(add$Timepoint == "BL Redo"), 2] = "Pre-treatment"
add[which(add$Timepoint == "W2"), 2] = "2 weeks"
add[which(add$Timepoint == "w2"), 2] = "2 weeks"
add[which(add$Timepoint == "w2 Redo"), 2] = "2 weeks"
add[which(add$Timepoint == "M"), 2] = "4 weeks/Mid-treatment"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[add$Timepoint == "Pre-treatment"] = "T1"
add$Timepoint_coded[add$Timepoint == "2 weeks"] = "T2"
add$Timepoint_coded[add$Timepoint == "4 weeks/Mid-treatment"] = "T2.5"
Entrez_pheno[[8]] = Entrez_pheno[[8]] %>% inner_join(add, by = "Sample.ID")

rm(add, pdata8)

# E4_1
# The set includes additional samples from the POETIC study and the NIT-2 study. 
# Here we will remove the samples from the NIT-2 study that have been reanalysed,
# because these are no-intervention treatment samples from GSE73237. The samples that
# have later been reanalysed in GSE126870 will also be removed.

# Since response is not reported, all of the samples that we are keeping here
# are gonna have missing values on response.

pdata9 = pData(GEOsets[["E4_1"]])
add = pdata9 %>%
  dplyr::select(`geo_accession`, relation, description, `group:ch1`, `timepoint:ch1`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Relation = relation, Desc = description, 
                Treated = `group:ch1`, Timepoint_pre = `timepoint:ch1`) %>%
  dplyr::filter(!grepl("Reanalyzed by: ", Relation)) %>%
  dplyr::filter(!grepl("Reanalysis of: ", Relation)) %>%
  dplyr::select(-Relation, -Treated) # Treated is removed because all subjects were treated
add$Timepoint = add$Timepoint_pre
add$Timepoint[add$Timepoint_pre == "baseline (B)"] = "Pre-treatment"
add$Timepoint[add$Timepoint_pre == "surgery (S)"] = "2 weeks (surgery)"
add$Timepoint_coded = add$Timepoint
add$Timepoint_coded[add$Timepoint == "Pre-treatment"] = "T1"
add$Timepoint_coded[add$Timepoint == "2 weeks (surgery)"] = "T2"
add$Desc = gsub("Treated.", "", add$Desc)
add$Desc = gsub("B", "", add$Desc)
add$Desc = gsub("S", "", add$Desc)
str_sub(add$Desc, 0, 0) = "E4."
add = add %>% dplyr::select(-Timepoint_pre) %>%
  dplyr::rename(Patient.ID = Desc)

# E4_2
# We are going to use the official supplementary file too here, because unlike
# the pheno info deposited in GEO, it contains info on response.
# Only baseline samples (reanalysed or not) are included here. These will be
# matched to the surgery samples from E4_1. However, the info about which study
# each sample came from will be retained to be included as a batch effect

pdata10 = pData(GEOsets[["E4_2"]])
pdata10_ext = openxlsx::read.xlsx("data/Supplementary/E4_2_supp.xlsx", sheet = 5,
                                  startRow = 3)
add1 = pdata10 %>%
  dplyr::select(`geo_accession`, title, `timepoint:ch1`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Patient.ID = title, Timepoint = `timepoint:ch1`)
add1$Patient.ID = gsub("AI.Treated.", "", add1$Patient.ID)
add1$Patient.ID = gsub(" \\[reanalysis\\]", "", add1$Patient.ID)
add1$Patient.ID = gsub("B", "", add1$Patient.ID)
str_sub(add1$Patient.ID, 0, 0) = "E4."
add1$Timepoint = gsub("baseline \\(B\\)", "Pre-treatment", add1$Timepoint)
add1$Timepoint_coded = "T1"

pdata10_ext = pdata10_ext[, c(2,11)]
colnames(pdata10_ext) = c("Patient.ID", "Response")
add2 = pdata10_ext %>% # Patient and Response columns
  dplyr::filter(!grepl("Control.", Patient.ID))
add2$Response = gsub("responder", "Responder", add2$Response)
add2$Response = gsub("non-Responder", "Non_responder", add2$Response)
add2$Response_type = "Ki67 60% expr. decr."
add2$Patient.ID = gsub("Treated", "E4", add2$Patient.ID)

add3 = add1 %>% left_join(add2, by = "Patient.ID") # for E4_2
add4 = add %>% left_join(add2, by = "Patient.ID") # for E4_1

Entrez_pheno[[9]] = Entrez_pheno[[9]] %>% inner_join(add4, by = "Sample.ID")
Entrez_pheno[[10]] = Entrez_pheno[[10]] %>% inner_join(add3, by = "Sample.ID")

rm(add, add1, add2, add3, add4, pdata10, pdata10_ext, pdata9)

# Adding a coded response column in each dataset
# C1 - using a supplementary data file
C1_Response = openxlsx::read.xlsx("data/Supplementary/C1.Response.xlsx", sheet = 2, colNames = T)
Entrez_pheno[[1]] = Entrez_pheno[[1]] %>% inner_join(C1_Response, by = "Patient.ID")
Entrez_pheno[[1]]$Response[Entrez_pheno[[1]]$Response_coded == 1] = "Responder"
Entrez_pheno[[1]]$Response[Entrez_pheno[[1]]$Response_coded == 0] = "Non_responder"

rm(C1_Response)

# C2
Entrez_pheno[[2]]$Response_coded = NA
Entrez_pheno[[2]]$Response_coded[Entrez_pheno[[2]]$Response == "Responder"] = 1
Entrez_pheno[[2]]$Response_coded[Entrez_pheno[[2]]$Response == "Non_responder"] = 0

# C3
Entrez_pheno[[3]]$Response_coded = NA
Entrez_pheno[[3]]$Response_coded[Entrez_pheno[[3]]$Response == "Responder"] = 1
Entrez_pheno[[3]]$Response_coded[Entrez_pheno[[3]]$Response == "Non_responder"] = 0

# E1_1
Entrez_pheno[[4]]$Response_coded = NA
Entrez_pheno[[4]]$Response_coded[Entrez_pheno[[4]]$Response == "Responder"] = 1
Entrez_pheno[[4]]$Response_coded[Entrez_pheno[[4]]$Response == "Non_responder"] = 0

# E1_2 (all are responders)
Entrez_pheno[[5]]$Response_coded = 1

# E1_3
Entrez_pheno[[6]]$Response_coded = NA
Entrez_pheno[[6]]$Response_coded[Entrez_pheno[[6]]$Response == "Responder"] = 1
Entrez_pheno[[6]]$Response_coded[Entrez_pheno[[6]]$Response == "Non_responder"] = 0

# E2
Entrez_pheno[[7]]$Response_coded = NA
Entrez_pheno[[7]]$Response_coded[Entrez_pheno[[7]]$Response == "Responder"] = 1
Entrez_pheno[[7]]$Response_coded[Entrez_pheno[[7]]$Response == "Non_responder"] = 0

# E3
Entrez_pheno[[8]]$Response_coded = NA
Entrez_pheno[[8]]$Response_coded[Entrez_pheno[[8]]$Response == "Responder"] = 1
Entrez_pheno[[8]]$Response_coded[Entrez_pheno[[8]]$Response == "Non_responder"] = 0

# E4_1
Entrez_pheno[[9]]$Response_coded = NA
Entrez_pheno[[9]]$Response_coded[Entrez_pheno[[9]]$Response == "Responder"] = 1
Entrez_pheno[[9]]$Response_coded[Entrez_pheno[[9]]$Response == "Non_responder"] = 0

# E4_2
Entrez_pheno[[10]]$Response_coded = NA
Entrez_pheno[[10]]$Response_coded[Entrez_pheno[[10]]$Response == "Responder"] = 1
Entrez_pheno[[10]]$Response_coded[Entrez_pheno[[10]]$Response == "Non_responder"] = 0

# Joining in on large phenotypic dataset
for (i in 1:length(Entrez_pheno)){
  Entrez_pheno[[i]]$Treatment_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "No_intervention", 0, 1)
  Entrez_pheno[[i]]$Endo_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "Endocrine_treatment", 1, 0)
  Entrez_pheno[[i]]$Chemo_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "Chemotherapy", 1, 0)
  Entrez_pheno[[i]] = Entrez_pheno[[i]] %>% 
    dplyr::select(Sample.ID, Patient.ID, Dataset, Treatment, Treatment_status, 
                  Chemo_status, Endo_status, Response_type, Response, 
                  Response_coded, Timepoint, Timepoint_coded)
}

Pheno = rbind(Entrez_pheno[[1]], Entrez_pheno[[2]], Entrez_pheno[[3]],
              Entrez_pheno[[4]], Entrez_pheno[[5]], Entrez_pheno[[6]],
              Entrez_pheno[[7]], Entrez_pheno[[8]], Entrez_pheno[[9]],
              Entrez_pheno[[10]]); rm(i)

# Platform metadata
Pheno$Platform = as.character(Pheno$Dataset)
Pheno$Platform[Pheno$Dataset == "C1"] = "GPL17303"
Pheno$Platform[Pheno$Dataset == "C2"] = "GPL14668"
Pheno$Platform[Pheno$Dataset == "C3"] = "GPL16791"
Pheno$Platform[Pheno$Dataset == "E1_1"] = "GPL96"
Pheno$Platform[Pheno$Dataset == "E1_2"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E1_3"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E2"] = "GPL13376"
Pheno$Platform[Pheno$Dataset == "E3"] = "GPL6480"
Pheno$Platform[Pheno$Dataset == "E4_1"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E4_2"] = "GPL10558"

Pheno$Platform_comp = as.character(Pheno$Platform)
Pheno$Platform_comp[Pheno$Dataset == "C1"] = "Ampliseq"
Pheno$Platform_comp[Pheno$Dataset == "C2"] = "Agilent"
Pheno$Platform_comp[Pheno$Dataset == "C3"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E1_1"] = "Affymetrix"
Pheno$Platform_comp[Pheno$Dataset == "E1_2"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E1_3"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E2"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E3"] = "Agilent"
Pheno$Platform_comp[Pheno$Dataset == "E4_1"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E4_2"] = "Illumina"

rm(Entrez_pheno); gc()

Pheno$Response = factor(Pheno$Response, levels = c("Responder", "Non_responder"),
                        labels = c("Responder", "Non_responder"))
Pheno$Treatment = factor(Pheno$Treatment, levels = c("Chemotherapy", "Endocrine_treatment"),
                         labels = c("Chemotherapy", "Endocrine_treatment"))
Pheno$Timepoint_coded = factor(Pheno$Timepoint_coded, levels = c("T1", "T2", "T2.5", "T3"),
                               labels = c("T1", "T2", "T2.5", "T3"))
Pheno$Dataset = factor(Pheno$Dataset, levels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                                                 "E4_1", "E4_2"),
                       labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                                  "E4_1", "E4_2"))
Pheno_response = Pheno %>% dplyr::filter(is.na(Response) == FALSE)
Pheno_times = Pheno %>% dplyr::filter(Timepoint_coded == "T1" | Timepoint_coded == "T2")
Pheno_exprs = Pheno_response %>% 
  dplyr::filter(Timepoint_coded == "T1" | Timepoint_coded == "T2")

# Stratified splitting into training, validation and test #####
RNGversion("4.0.2")
set.seed(123)
tobesplit = Pheno_exprs %>%
  mutate(nr = row_number()) %>%
  dplyr::select(nr, everything()) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
train_set = tobesplit %>%
  group_by(Response, Timepoint_coded, Dataset) %>%
  sample_frac(0.7) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Response, Timepoint_coded, Dataset) %>%
  sample_frac(0.65) %>%
  as.data.frame()
test_set = anti_join(tobesplit, as.data.frame(rbind(train_set, validation_set)))

train_set = train_set %>% dplyr::ungroup()
validation_set = validation_set %>% dplyr::ungroup()
test_set = test_set %>% dplyr::ungroup()

train_samples = train_set$Sample.ID
validation_samples = validation_set$Sample.ID
test_samples = test_set$Sample.ID
trainval_samples = c(train_samples, validation_samples)

GEOsets_trainval = list(); GEOsets_test = list()
# For the test sets in the first loop we also add the first sample of each GEOset
# object which was dropped in lines 60-68
for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_trainval[[i]] = GEOsets[[i]][, intersect(trainval_samples, colnames(GEOsets[[i]]))]
  GEOsets_test[[i]] = GEOsets[[i]][, c(colnames(GEOsets[[i]])[1],
                                       intersect(test_samples, colnames(GEOsets[[i]])))]
}
for (i in c(1, 3, 7)){
  GEOsets_trainval[[i]] = GEOsets[[i]][, c("X", 
                                           intersect(trainval_samples, colnames(GEOsets[[i]])))]
  GEOsets_test[[i]] = GEOsets[[i]][, c("X", 
                                       intersect(test_samples, colnames(GEOsets[[i]])))]
}
names(GEOsets_test) = names(GEOsets_trainval) = names(GEOsets)
rm(i); gc()

# Probe removal (missing values) and kNN-imputation #####
# Rows that contain at least 25% of missing values are filtered out

for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_trainval[[i]] = GEOsets_trainval[[i]][rowSums(is.na(GEOsets_trainval[[i]]@assayData[["exprs"]]))/
                                length(colnames(GEOsets_trainval[[i]]@assayData[["exprs"]])) < 0.25, ]
}
rm(i)
for(i in c(1, 3, 7)){
  matrix = as.matrix(GEOsets_trainval[[i]][,2:ncol(GEOsets_trainval[[i]])])
  rownames(matrix) = GEOsets_trainval[[i]]$X
  class(matrix) = "numeric"
  matrix = matrix[rowSums(is.na(matrix))/
                    length(colnames(matrix)) < 0.25, ]
  GEOsets_trainval[[i]] = as.data.frame(matrix)
  GEOsets_trainval[[i]]$X = rownames(matrix)
  GEOsets_trainval[[i]] = GEOsets_trainval[[i]] %>% dplyr::select(X, everything())
}
rm(matrix)

# The rest of the missing values are imputed using kNN for k = 10
GEO_knn_pre1 = list(GEOsets_trainval[["C2"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_1"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_2"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_3"]]@assayData[["exprs"]],
               GEOsets_trainval[["E3"]]@assayData[["exprs"]],
               GEOsets_trainval[["E4_1"]]@assayData[["exprs"]],
               GEOsets_trainval[["E4_2"]]@assayData[["exprs"]])

GEO_knn_pre2 = list(GEOsets_trainval[["C1"]],
                    GEOsets_trainval[["C3"]],
                    GEOsets_trainval[["E2"]])

# kNN
for(i in 1:length(GEO_knn_pre1)){
  RNGversion("4.0.2")
  GEO_knn_pre1[[i]] = impute.knn(GEO_knn_pre1[[i]], k = 10, maxp = nrow(GEO_knn_pre1[[i]]),
                            rng.seed = 123)
  GEO_knn_pre1[[i]] = GEO_knn_pre1[[i]][["data"]]
}
rm(i)
for(i in 1:length(GEO_knn_pre2)){
  matrix = as.matrix(GEO_knn_pre2[[i]][,2:ncol(GEO_knn_pre2[[i]])])
  rownames(matrix) = GEO_knn_pre2[[i]]$X
  class(matrix) = "numeric"
  RNGversion("4.0.2")
  matrix = impute.knn(matrix, k = 10, maxp = nrow(matrix),
                      rng.seed = 123)[["data"]]
  GEO_knn_pre2[[i]] = as.data.frame(matrix)
  GEO_knn_pre2[[i]]$X = rownames(matrix)
  GEO_knn_pre2[[i]] = GEO_knn_pre2[[i]] %>% dplyr::select(X, everything())
}
rm(i, matrix)
names(GEO_knn_pre1) = names(GEOsets)[c(2, 4, 5, 6, 8, 9, 10)]
names(GEO_knn_pre2) = names(GEOsets)[c(1, 3, 7)]

# In order to make the test set consistent with the trainval sets we need to
# ensure the same probes are kept per study and appropriate kNN-imputation is
# performed:

for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_test[[i]] = GEOsets_test[[i]][rownames(GEOsets_trainval[[i]]), ]
}
rm(i)
for(i in c(1, 3, 7)){
  GEOsets_test[[i]] = GEOsets_test[[i]] %>%
    dplyr::filter(X %in% intersect(GEOsets_trainval[[i]]$X, GEOsets_test[[i]]$X))
}
rm(i)

# In order to impute the test samples for each study, we do the following in order:
# 1. find the 10 nearest neighbors for every gene in our imputed training and validation
#    samples of a study.
# 2. for each sample, impute the missing value in a gene using 10-NN averaging of the 
#    values in that gene's neighbours, as these were determined in step 1.
# 3. repeat for all studies

# For the external validation sets, we normalise within each study later

# The following processes are very expensive memory-wise and require a supercluster
# to run. Locally, they can be run one at a time. The loop needs excessive memory

neighbors_list_1 = list()
for (i in 1:length(GEO_knn_pre1)){
  neighbors = as.data.frame(matrix(0, ncol = 1, nrow = 10))
  colnames(neighbors) = "scaffold"
  dists = dist(GEO_knn_pre1[[i]], method = "euclidean"); gc()
  dists2 = as.matrix(dists); rm(dists); gc()
  for (j in 1:nrow(dists2)){
    gene_j_neighbors = names(sort(dists2[j,], decreasing = FALSE)[2:11]) # 2:11 is used for the top 10 neighbor genes, because the first one will always be the gene itself
    neighbors = cbind(neighbors, gene_j_neighbors)
  }
  gc()
  neighbors = neighbors[,-1]
  colnames(neighbors) = rownames(dists2)
  neighbors_list_1[[i]] = neighbors
  rm(neighbors, dists2); gc()
}
names(neighbors_list_1) = names(GEO_knn_pre1)

neighbors_list_2 = list()
for (i in 1:length(GEO_knn_pre2)){
  neighbors = as.data.frame(matrix(0, ncol = 1, nrow = 10))
  colnames(neighbors) = "scaffold"
  gem = as.matrix(GEO_knn_pre2[[i]][,-1]); rownames(gem) = GEO_knn_pre2[[i]]$X
  dists = dist(gem, method = "euclidean"); gc()
  dists2 = as.matrix(dists); rm(dists); gc()
  for (j in 1:nrow(dists2)){
    gene_j_neighbors = names(sort(dists2[j,], decreasing = FALSE)[2:11]) # 2:11 is used for the top 10 neighbor genes, because the first one will always be the gene itself
    neighbors = cbind(neighbors, gene_j_neighbors)
  }
  gc()
  neighbors = neighbors[,-1]
  colnames(neighbors) = rownames(dists2)
  neighbors_list_2[[i]] = neighbors
  rm(neighbors, dists2); gc()
}
names(neighbors_list_2) = names(GEO_knn_pre2)

# Imputing the sets that correspond to neighbors_list_1:
GEOsets_test_new = list()
for (i in c(2, 4, 5, 6, 8, 9, 10)) {
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = rownames(GEOsets_test[[i]])
  imp_map = neighbors_list_1[[names(GEOsets_test)[i]]]
  for (j in 1:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]]@assayData[["exprs"]][, j]
    for (k in 1:length(sample)) {
      if (is.na(sample[k]) == TRUE) {
        probe = names(sample[k])
        neighbors = imp_map[, probe]
        neighbor_values = sample[neighbors]
        if (sum(is.na(neighbor_values)) == 10) {
          sample[k] = mean(sample, na.rm = TRUE)
        } else {
          sample[k] = mean(neighbor_values, na.rm = TRUE)
        }
      }
    }
    scaffold = cbind(scaffold, sample)
    rm(sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  class(GEOsets_test_new[[i]]) = "numeric"
  rownames(GEOsets_test_new[[i]]) = rownames(GEOsets_test[[i]])
  rm(scaffold, probe, neighbor_values, neighbors)
}
rm(i, j, k); gc()

# Imputing the sets that correspond to neighbors_list_2:
for (i in c(1, 3, 7)) {
  scaffold = matrix(GEOsets_test[[i]]$X)
  colnames(scaffold) = "X"
  imp_map = neighbors_list_2[[names(GEOsets_test)[i]]]
  for (j in 2:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]][, j]
    names(sample) = GEOsets_test[[i]]$X
    for (k in 1:length(sample)) {
      if (is.na(sample[k]) == TRUE) {
        probe = names(sample[k])
        neighbors = imp_map[, probe]
        neighbor_values = sample[neighbors]
        if (sum(is.na(neighbor_values)) == 10) {
          sample[k] = mean(sample, na.rm = TRUE)
        } else {
          sample[k] = mean(neighbor_values, na.rm = TRUE)
        }
      }
    }
    scaffold = cbind(scaffold, sample)
    rm(sample)
  }
  GEOsets_test_new[[i]] = as.data.frame(scaffold)
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  rownames(GEOsets_test_new[[i]]) = GEOsets_test_new[[i]]$X
  GEOsets_test_new[[i]][colnames(GEOsets_test_new[[i]])[2:ncol(GEOsets_test_new[[i]])]] = 
   lapply(GEOsets_test_new[[i]][colnames(GEOsets_test_new[[i]])[2:ncol(GEOsets_test_new[[i]])]], 
         function(x) as.numeric(x))
  rm(scaffold, probe, neighbor_values, neighbors)
}
rm(i, j, k, imp_map, gene_j_neighbors); gc()
names(GEOsets_test_new) = names(GEOsets_test)

################################################################################


# The following code can be used if the test samples are to be imputed one 
# by one through an "add a sample -> impute -> remove sample" manner

# Second test set imputation approach - NOT USED IN OUR CASE #####
# Imputing the test samples one by one (for each study separately)
GEOsets_test_new = list()
for (i in c(2, 4, 5, 6, 8, 9, 10)){
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = rownames(GEOsets_test[[i]])
  impset = GEO_knn_pre1[[names(GEOsets_test)[i]]]
  for (j in 1:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]]@assayData[["exprs"]][, j]
    impmatrix = cbind(impset, sample)
    RNGversion("4.0.2")
    imp = impute.knn(impmatrix, k = 10, maxp = nrow(impmatrix),
                                   rng.seed = 123)
    imp = imp[["data"]]
    imp_sample = imp[, ncol(imp)]
    scaffold = cbind(scaffold, imp_sample)
    rm(sample, impmatrix, imp, imp_sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  class(GEOsets_test_new[[i]]) = "numeric"
  rm(scaffold, impset)
}
rm(i, j)

for (i in c(1, 3, 7)){
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = GEOsets_test[[i]]$X
  impset = as.matrix(GEO_knn_pre2[[names(GEOsets_test)[i]]][,2:ncol(GEO_knn_pre2[[names(GEOsets_test)[i]]])])
  rownames(impset) = GEO_knn_pre2[[names(GEOsets_test)[i]]]$X
  class(impset) = "numeric"
  for (j in 1:(ncol(GEOsets_test[[i]]) - 1)) {
    sample = GEOsets_test[[i]][, (j + 1)]
    impmatrix = cbind(impset, sample)
    RNGversion("4.0.2")
    imp = impute.knn(impmatrix, k = 10, maxp = nrow(impmatrix),
                     rng.seed = 123)
    imp = imp[["data"]]
    imp_sample = imp[, ncol(imp)]
    scaffold = cbind(scaffold, imp_sample)
    rm(sample, impmatrix, imp, imp_sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  X = rownames(GEOsets_test_new[[i]])
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])[2:ncol(GEOsets_test[[i]])]
  GEOsets_test_new[[i]] = as.data.frame(GEOsets_test_new[[i]])
  GEOsets_test_new[[i]]$X = X; rm(X)
  GEOsets_test_new[[i]] = GEOsets_test_new[[i]] %>% dplyr::select(X, everything())
  rm(scaffold, impset)
}

rm(i, j)
names(GEOsets_test_new) = names(GEOsets_test)

################################################################################

# Annotation: files and function ######
# ID_Map
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]

alias = org.Hs.egALIAS2EG
mapped_genes_alias = mappedkeys(alias)
alias_df = as.data.frame(alias[mapped_genes_alias])
alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
alias_df$HGNC_Official = "No"

ID_Map = rbind(official_df, alias_df) %>% distinct()
ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
  dplyr::rename(probe=Gene.Symbol) %>%
  dplyr::select(probe, EntrezGene.ID, HGNC_Official)

# Aliases
aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
Aliases = official_df %>% inner_join(aliases_for_join,
                                     by = "EntrezGene.ID") %>%
  dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
  dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
                Entrez = EntrezGene.ID) %>%
  distinct()

rm(alias, alias_df, aliases_for_join, official, mapped_genes_alias, mapped_genes_official)

# Load the annotation and genefu predictions functions
source("R_scripts/Genefu_functions.R")

# Annotation: loaded sets #####
# Annotating using annot() on the imputed data
Loaded_sets_trainval = list(GEO_knn_pre2[["C1"]], GEO_knn_pre2[["C3"]],
                            GEO_knn_pre2[["E2"]])
names(Loaded_sets_trainval) = c("C1", "C3", "E2")
Loaded_sets_test = list(GEOsets_test_new[["C1"]], GEOsets_test_new[["C3"]],
                   GEOsets_test_new[["E2"]])
names(Loaded_sets_test) = c("C1", "C3", "E2")

Loaded_exprs_Entrez_trainval = list()
Loaded_exprs_HGNC_trainval = list()
Entrez_dannot_trainval = list()
HGNC_dannot_trainval = list()

Loaded_exprs_Entrez_test = list()
Loaded_exprs_HGNC_test = list()
Entrez_dannot_test = list()
HGNC_dannot_test = list()

for (i in 1:length(Loaded_sets_trainval)){
  # Entrez-annotated data frames
  annotation = annot(expressionDF = Loaded_sets_trainval[[i]], 
                     ID_Map = ID_Map, Aliases = Aliases)
  Loaded_exprs_Entrez_trainval[[i]] = annotation$Entrez
  Loaded_exprs_Entrez_trainval[[i]] = Loaded_exprs_Entrez_trainval[[i]] %>% 
    dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
    group_by(EntrezGene.ID) %>%
    summarise_all(mean, na.rm = TRUE)
  
  Entrez_dannot_trainval[[i]] = as.data.frame(Loaded_exprs_Entrez_trainval[[i]][,"EntrezGene.ID"])
  Entrez_dannot_trainval[[i]]$probe = Entrez_dannot_trainval[[i]]$EntrezGene.ID
  Entrez_dannot_trainval[[i]] = Entrez_dannot_trainval[[i]] %>%
    dplyr::select(probe, EntrezGene.ID)
  rownames(Entrez_dannot_trainval[[i]]) = Entrez_dannot_trainval[[i]]$probe
  
  # HGNC-annotated data frames
  Loaded_exprs_HGNC_trainval[[i]] = Loaded_exprs_Entrez_trainval[[i]] %>% inner_join(official_df, by = "EntrezGene.ID") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
    dplyr::filter(is.na(Gene.Symbol) == FALSE) %>%
    dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
    group_by(Gene.Symbol) %>%
    summarise_all(mean, na.rm = TRUE)
  
  HGNC_dannot_trainval[[i]] = as.data.frame(Loaded_exprs_HGNC_trainval[[i]][,"Gene.Symbol"])
  HGNC_dannot_trainval[[i]]$probe = HGNC_dannot_trainval[[i]]$Gene.Symbol
  HGNC_dannot_trainval[[i]] = HGNC_dannot_trainval[[i]] %>%
    dplyr::select(probe, Gene.Symbol)
  rownames(HGNC_dannot_trainval[[i]]) = HGNC_dannot_trainval[[i]]$probe
}

for (i in 1:length(Loaded_sets_test)){
  # Entrez-annotated data frames
  annotation = annot(expressionDF = Loaded_sets_test[[i]], 
                     ID_Map = ID_Map, Aliases = Aliases)
  Loaded_exprs_Entrez_test[[i]] = annotation$Entrez
  Loaded_exprs_Entrez_test[[i]] = Loaded_exprs_Entrez_test[[i]] %>% 
    dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
    group_by(EntrezGene.ID) %>%
    summarise_all(mean, na.rm = TRUE)
  
  Entrez_dannot_test[[i]] = as.data.frame(Loaded_exprs_Entrez_test[[i]][,"EntrezGene.ID"])
  Entrez_dannot_test[[i]]$probe = Entrez_dannot_test[[i]]$EntrezGene.ID
  Entrez_dannot_test[[i]] = Entrez_dannot_test[[i]] %>%
    dplyr::select(probe, EntrezGene.ID)
  rownames(Entrez_dannot_test[[i]]) = Entrez_dannot_test[[i]]$probe
  
  # HGNC-annotated data frames
  Loaded_exprs_HGNC_test[[i]] = Loaded_exprs_Entrez_test[[i]] %>% inner_join(official_df, by = "EntrezGene.ID") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
    dplyr::filter(is.na(Gene.Symbol) == FALSE) %>%
    dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
    group_by(Gene.Symbol) %>%
    summarise_all(mean, na.rm = TRUE)
  
  HGNC_dannot_test[[i]] = as.data.frame(Loaded_exprs_HGNC_test[[i]][,"Gene.Symbol"])
  HGNC_dannot_test[[i]]$probe = HGNC_dannot_test[[i]]$Gene.Symbol
  HGNC_dannot_test[[i]] = HGNC_dannot_test[[i]] %>%
    dplyr::select(probe, Gene.Symbol)
  rownames(HGNC_dannot_test[[i]]) = HGNC_dannot_test[[i]]$probe
}
rm(i, pam50, annotation); gc()
names(Loaded_exprs_Entrez_test) = names(Loaded_exprs_Entrez_trainval) = 
  names(Loaded_exprs_HGNC_test) = names(Loaded_exprs_HGNC_trainval) = 
  names(Entrez_dannot_test) = names(Entrez_dannot_trainval) = 
  names(HGNC_dannot_test) = names(HGNC_dannot_trainval) = 
  names(Loaded_sets_test)

# Annotation: sets from GEOquery (sets with whole experiments in the environment) #####
# 5 kinds of objects are created here for each set: 
# 1) exprs# is the clean expression matrix (probes matching to unique genes) - deleted later, 
# 2) exprs#_Entrez is the clean expression matrix with Entrez ID's
# 3) dannot#_Entrez is the Entrez annotation data frame for that matrix,
# 4) exprs#_HGNC is the clean expression matrix with gene symbols and
# 5) dannot#_HGNC is the HGNC annotation data frame for that second matrix

# GSE32603 - C2
exprs32603_trainval = as.data.frame(GEO_knn_pre1[["C2"]])
exprs32603_trainval$X = GEOsets_trainval[["C2"]]@featureData@data[["GENE SYMBOL"]]
exprs32603_trainval = exprs32603_trainval %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID
exprs32603_Entrez_trainval = annot(expressionDF = exprs32603_trainval, ID_Map = ID_Map,
                          Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID

dannot32603_Entrez_trainval = as.data.frame(exprs32603_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot32603_Entrez_trainval) = "EntrezGene.ID"
dannot32603_Entrez_trainval$probe = dannot32603_Entrez_trainval$EntrezGene.ID
dannot32603_Entrez_trainval = dannot32603_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot32603_Entrez_trainval) = dannot32603_Entrez_trainval$probe

# In the following chunck of code we use official_df to map Entrez ID's to official 
# Gene Symbols which are used by the Integrative Cluster algorithm. Aliases could mislead
# the algorithm in terms of ID's

exprs32603_HGNC_trainval = exprs32603_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot32603_HGNC_trainval = as.data.frame(exprs32603_HGNC_trainval[,"Gene.Symbol"])
dannot32603_HGNC_trainval$probe = dannot32603_HGNC_trainval$Gene.Symbol
dannot32603_HGNC_trainval = dannot32603_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot32603_HGNC_trainval) = dannot32603_HGNC_trainval$probe

rm(exprs32603_trainval)

# test set samples
exprs32603_test = as.data.frame(GEOsets_test_new[["C2"]])
exprs32603_test$X = GEOsets_test[["C2"]]@featureData@data[["GENE SYMBOL"]]
exprs32603_test = exprs32603_test %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID
exprs32603_Entrez_test = annot(expressionDF = exprs32603_test, ID_Map = ID_Map,
                               Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID

dannot32603_Entrez_test = as.data.frame(exprs32603_Entrez_test[,"EntrezGene.ID"])
colnames(dannot32603_Entrez_test) = "EntrezGene.ID"
dannot32603_Entrez_test$probe = dannot32603_Entrez_test$EntrezGene.ID
dannot32603_Entrez_test = dannot32603_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot32603_Entrez_test) = dannot32603_Entrez_test$probe

# In the following chunck of code we use official_df to map Entrez ID's to official 
# Gene Symbols which are used by the Integrative Cluster algorithm. Aliases could mislead
# the algorithm in terms of ID's

exprs32603_HGNC_test = exprs32603_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot32603_HGNC_test = as.data.frame(exprs32603_HGNC_test[,"Gene.Symbol"])
dannot32603_HGNC_test$probe = dannot32603_HGNC_test$Gene.Symbol
dannot32603_HGNC_test = dannot32603_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot32603_HGNC_test) = dannot32603_HGNC_test$probe

rm(exprs32603_test)

# GSE20181 - E1_1
feature_trainval = GEOsets_trainval[["E1_1"]]@featureData@data %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
exprs20181_trainval = as.data.frame(GEO_knn_pre1[["E1_1"]])[feature_trainval$ID,]
exprs20181_trainval$ID = rownames(exprs20181_trainval)
exprs20181_Entrez_trainval = exprs20181_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  dplyr::filter(is.na(ENTREZ_GENE_ID) == FALSE) %>% # just in case
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = ENTREZ_GENE_ID)

dannot20181_Entrez_trainval = as.data.frame(exprs20181_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot20181_Entrez_trainval) = "EntrezGene.ID"
dannot20181_Entrez_trainval$probe = dannot20181_Entrez_trainval$EntrezGene.ID
dannot20181_Entrez_trainval = dannot20181_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot20181_Entrez_trainval) = dannot20181_Entrez_trainval$probe

exprs20181_HGNC_trainval = exprs20181_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot20181_HGNC_trainval = as.data.frame(exprs20181_HGNC_trainval[,"Gene.Symbol"])
dannot20181_HGNC_trainval$probe = dannot20181_HGNC_trainval$Gene.Symbol
dannot20181_HGNC_trainval = dannot20181_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot20181_HGNC_trainval) = dannot20181_HGNC_trainval$probe

rm(exprs20181_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_1"]]@featureData@data %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
exprs20181_test = as.data.frame(GEOsets_test_new[["E1_1"]])[feature_test$ID,]
exprs20181_test$ID = rownames(exprs20181_test)
exprs20181_Entrez_test = exprs20181_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  dplyr::filter(is.na(ENTREZ_GENE_ID) == FALSE) %>% # just in case
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = ENTREZ_GENE_ID)

dannot20181_Entrez_test = as.data.frame(exprs20181_Entrez_test[,"EntrezGene.ID"])
colnames(dannot20181_Entrez_test) = "EntrezGene.ID"
dannot20181_Entrez_test$probe = dannot20181_Entrez_test$EntrezGene.ID
dannot20181_Entrez_test = dannot20181_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot20181_Entrez_test) = dannot20181_Entrez_test$probe

exprs20181_HGNC_test = exprs20181_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot20181_HGNC_test = as.data.frame(exprs20181_HGNC_test[,"Gene.Symbol"])
dannot20181_HGNC_test$probe = dannot20181_HGNC_test$Gene.Symbol
dannot20181_HGNC_test = dannot20181_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot20181_HGNC_test) = dannot20181_HGNC_test$probe

rm(exprs20181_test, feature_test)

# GSE55374 - E1_2
feature_trainval = GEOsets_trainval[["E1_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs55374_trainval = as.data.frame(GEO_knn_pre1[["E1_2"]])[feature_trainval$ID,]
exprs55374_trainval$ID = rownames(exprs55374_trainval)
exprs55374_Entrez_trainval = exprs55374_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs55374_Entrez_trainval$EntrezGene.ID = as.character(exprs55374_Entrez_trainval$EntrezGene.ID)

dannot55374_Entrez_trainval = as.data.frame(exprs55374_Entrez_trainval[,"EntrezGene.ID"])
dannot55374_Entrez_trainval$probe = dannot55374_Entrez_trainval$EntrezGene.ID
dannot55374_Entrez_trainval = dannot55374_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot55374_Entrez_trainval) = dannot55374_Entrez_trainval$probe

exprs55374_HGNC_trainval = exprs55374_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot55374_HGNC_trainval = as.data.frame(exprs55374_HGNC_trainval[,"Gene.Symbol"])
dannot55374_HGNC_trainval$probe = dannot55374_HGNC_trainval$Gene.Symbol
dannot55374_HGNC_trainval = dannot55374_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot55374_HGNC_trainval) = dannot55374_HGNC_trainval$probe

rm(exprs55374_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs55374_test = as.data.frame(GEOsets_test_new[["E1_2"]])[feature_test$ID,]
exprs55374_test$ID = rownames(exprs55374_test)
exprs55374_Entrez_test = exprs55374_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs55374_Entrez_test$EntrezGene.ID = as.character(exprs55374_Entrez_test$EntrezGene.ID)

dannot55374_Entrez_test = as.data.frame(exprs55374_Entrez_test[,"EntrezGene.ID"])
dannot55374_Entrez_test$probe = dannot55374_Entrez_test$EntrezGene.ID
dannot55374_Entrez_test = dannot55374_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot55374_Entrez_test) = dannot55374_Entrez_test$probe

exprs55374_HGNC_test = exprs55374_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot55374_HGNC_test = as.data.frame(exprs55374_HGNC_test[,"Gene.Symbol"])
dannot55374_HGNC_test$probe = dannot55374_HGNC_test$Gene.Symbol
dannot55374_HGNC_test = dannot55374_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot55374_HGNC_test) = dannot55374_HGNC_test$probe

rm(exprs55374_test, feature_test)

# GSE59515 - E1_3
feature_trainval = GEOsets_trainval[["E1_3"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs59515_trainval = as.data.frame(GEO_knn_pre1[["E1_3"]])[feature_trainval$ID,]
exprs59515_trainval$ID = rownames(exprs59515_trainval)
exprs59515_Entrez_trainval = exprs59515_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs59515_Entrez_trainval$EntrezGene.ID = as.character(exprs59515_Entrez_trainval$EntrezGene.ID)

dannot59515_Entrez_trainval = as.data.frame(exprs59515_Entrez_trainval[,"EntrezGene.ID"])
dannot59515_Entrez_trainval$probe = dannot59515_Entrez_trainval$EntrezGene.ID
dannot59515_Entrez_trainval = dannot59515_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot59515_Entrez_trainval) = dannot59515_Entrez_trainval$probe

exprs59515_HGNC_trainval = exprs59515_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot59515_HGNC_trainval = as.data.frame(exprs59515_HGNC_trainval[,"Gene.Symbol"])
dannot59515_HGNC_trainval$probe = dannot59515_HGNC_trainval$Gene.Symbol
dannot59515_HGNC_trainval = dannot59515_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot59515_HGNC_trainval) = dannot59515_HGNC_trainval$probe

rm(exprs59515_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_3"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs59515_test = as.data.frame(GEOsets_test_new[["E1_3"]])[feature_test$ID,]
exprs59515_test$ID = rownames(exprs59515_test)
exprs59515_Entrez_test = exprs59515_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs59515_Entrez_test$EntrezGene.ID = as.character(exprs59515_Entrez_test$EntrezGene.ID)

dannot59515_Entrez_test = as.data.frame(exprs59515_Entrez_test[,"EntrezGene.ID"])
dannot59515_Entrez_test$probe = dannot59515_Entrez_test$EntrezGene.ID
dannot59515_Entrez_test = dannot59515_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot59515_Entrez_test) = dannot59515_Entrez_test$probe

exprs59515_HGNC_test = exprs59515_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot59515_HGNC_test = as.data.frame(exprs59515_HGNC_test[,"Gene.Symbol"])
dannot59515_HGNC_test$probe = dannot59515_HGNC_test$Gene.Symbol
dannot59515_HGNC_test = dannot59515_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot59515_HGNC_test) = dannot59515_HGNC_test$probe

rm(exprs59515_test, feature_test)

# GSE87411 - E3
exprs87411_trainval = as.data.frame(GEO_knn_pre1[["E3"]])
exprs87411_trainval$X = GEOsets_trainval[["E3"]]@featureData@data[["GENE_SYMBOL"]]
exprs87411_trainval = exprs87411_trainval %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE)
exprs87411_Entrez_trainval = annot(expressionDF = exprs87411_trainval, ID_Map = ID_Map,
                          Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_Entrez_trainval = as.data.frame(exprs87411_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot87411_Entrez_trainval) = "EntrezGene.ID"
dannot87411_Entrez_trainval$probe = dannot87411_Entrez_trainval$EntrezGene.ID
dannot87411_Entrez_trainval = dannot87411_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot87411_Entrez_trainval) = dannot87411_Entrez_trainval$probe

exprs87411_HGNC_trainval = exprs87411_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_HGNC_trainval = as.data.frame(exprs87411_HGNC_trainval[,"Gene.Symbol"])
dannot87411_HGNC_trainval$probe = dannot87411_HGNC_trainval$Gene.Symbol
dannot87411_HGNC_trainval = dannot87411_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot87411_HGNC_trainval) = dannot87411_HGNC_trainval$probe

rm(exprs87411_trainval)

# test set samples
exprs87411_test = as.data.frame(GEOsets_test_new[["E3"]])
exprs87411_test$X = GEOsets_test[["E3"]]@featureData@data[["GENE_SYMBOL"]]
exprs87411_test = exprs87411_test %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE)
exprs87411_Entrez_test = annot(expressionDF = exprs87411_test, ID_Map = ID_Map,
                               Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_Entrez_test = as.data.frame(exprs87411_Entrez_test[,"EntrezGene.ID"])
colnames(dannot87411_Entrez_test) = "EntrezGene.ID"
dannot87411_Entrez_test$probe = dannot87411_Entrez_test$EntrezGene.ID
dannot87411_Entrez_test = dannot87411_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot87411_Entrez_test) = dannot87411_Entrez_test$probe

exprs87411_HGNC_test = exprs87411_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_HGNC_test = as.data.frame(exprs87411_HGNC_test[,"Gene.Symbol"])
dannot87411_HGNC_test$probe = dannot87411_HGNC_test$Gene.Symbol
dannot87411_HGNC_test = dannot87411_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot87411_HGNC_test) = dannot87411_HGNC_test$probe

rm(exprs87411_test)

# GSE105777 - E4_1
feature_trainval = GEOsets_trainval[["E4_1"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs105777_trainval = as.data.frame(GEO_knn_pre1[["E4_1"]])[feature_trainval$ID,]
exprs105777_trainval$ID = rownames(exprs105777_trainval)
exprs105777_Entrez_trainval = exprs105777_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs105777_Entrez_trainval$EntrezGene.ID = as.character(exprs105777_Entrez_trainval$EntrezGene.ID)

dannot105777_Entrez_trainval = as.data.frame(exprs105777_Entrez_trainval[,"EntrezGene.ID"])
dannot105777_Entrez_trainval$probe = dannot105777_Entrez_trainval$EntrezGene.ID
dannot105777_Entrez_trainval = dannot105777_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot105777_Entrez_trainval) = dannot105777_Entrez_trainval$probe

exprs105777_HGNC_trainval = exprs105777_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot105777_HGNC_trainval = as.data.frame(exprs105777_HGNC_trainval[,"Gene.Symbol"])
dannot105777_HGNC_trainval$probe = dannot105777_HGNC_trainval$Gene.Symbol
dannot105777_HGNC_trainval = dannot105777_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot105777_HGNC_trainval) = dannot105777_HGNC_trainval$probe

rm(exprs105777_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E4_1"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs105777_test = as.data.frame(GEOsets_test_new[["E4_1"]])[feature_test$ID,]
exprs105777_test$ID = rownames(exprs105777_test)
exprs105777_Entrez_test = exprs105777_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs105777_Entrez_test$EntrezGene.ID = as.character(exprs105777_Entrez_test$EntrezGene.ID)

dannot105777_Entrez_test = as.data.frame(exprs105777_Entrez_test[,"EntrezGene.ID"])
dannot105777_Entrez_test$probe = dannot105777_Entrez_test$EntrezGene.ID
dannot105777_Entrez_test = dannot105777_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot105777_Entrez_test) = dannot105777_Entrez_test$probe

exprs105777_HGNC_test = exprs105777_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot105777_HGNC_test = as.data.frame(exprs105777_HGNC_test[,"Gene.Symbol"])
dannot105777_HGNC_test$probe = dannot105777_HGNC_test$Gene.Symbol
dannot105777_HGNC_test = dannot105777_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot105777_HGNC_test) = dannot105777_HGNC_test$probe

rm(exprs105777_test, feature_test)

# GSE126870 - E4_2
feature_trainval = GEOsets_trainval[["E4_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs126870_trainval = as.data.frame(GEO_knn_pre1[["E4_2"]])[feature_trainval$ID,]
exprs126870_trainval$ID = rownames(exprs126870_trainval)
exprs126870_Entrez_trainval = exprs126870_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs126870_Entrez_trainval$EntrezGene.ID = as.character(exprs126870_Entrez_trainval$EntrezGene.ID)

dannot126870_Entrez_trainval = as.data.frame(exprs126870_Entrez_trainval[,"EntrezGene.ID"])
dannot126870_Entrez_trainval$probe = dannot126870_Entrez_trainval$EntrezGene.ID
dannot126870_Entrez_trainval = dannot126870_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot126870_Entrez_trainval) = dannot126870_Entrez_trainval$probe

exprs126870_HGNC_trainval = exprs126870_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot126870_HGNC_trainval = as.data.frame(exprs126870_HGNC_trainval[,"Gene.Symbol"])
dannot126870_HGNC_trainval$probe = dannot126870_HGNC_trainval$Gene.Symbol
dannot126870_HGNC_trainval = dannot126870_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot126870_HGNC_trainval) = dannot126870_HGNC_trainval$probe

rm(exprs126870_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E4_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs126870_test = as.data.frame(GEOsets_test_new[["E4_2"]])[feature_test$ID,]
exprs126870_test$ID = rownames(exprs126870_test)
exprs126870_Entrez_test = exprs126870_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs126870_Entrez_test$EntrezGene.ID = as.character(exprs126870_Entrez_test$EntrezGene.ID)

dannot126870_Entrez_test = as.data.frame(exprs126870_Entrez_test[,"EntrezGene.ID"])
dannot126870_Entrez_test$probe = dannot126870_Entrez_test$EntrezGene.ID
dannot126870_Entrez_test = dannot126870_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot126870_Entrez_test) = dannot126870_Entrez_test$probe

exprs126870_HGNC_test = exprs126870_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot126870_HGNC_test = as.data.frame(exprs126870_HGNC_test[,"Gene.Symbol"])
dannot126870_HGNC_test$probe = dannot126870_HGNC_test$Gene.Symbol
dannot126870_HGNC_test = dannot126870_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot126870_HGNC_test) = dannot126870_HGNC_test$probe

rm(exprs126870_test, feature_test)

# Joining them in a list of Entrez-annotated sets and a list of HGNC-annotated sets
Entrez_exprs_trainval = list(Loaded_exprs_Entrez_trainval[["C1"]], exprs32603_Entrez_trainval, 
                             Loaded_exprs_Entrez_trainval[["C3"]], exprs20181_Entrez_trainval, 
                             exprs55374_Entrez_trainval, exprs59515_Entrez_trainval, 
                    Loaded_exprs_Entrez_trainval[["E2"]], exprs87411_Entrez_trainval, 
                    exprs105777_Entrez_trainval, exprs126870_Entrez_trainval)
names(Entrez_exprs_trainval) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                        "E4_1", "E4_2")

Entrez_exprs_test = list(Loaded_exprs_Entrez_test[["C1"]], exprs32603_Entrez_test, 
                         Loaded_exprs_Entrez_test[["C3"]], exprs20181_Entrez_test, 
                         exprs55374_Entrez_test, exprs59515_Entrez_test, 
                         Loaded_exprs_Entrez_test[["E2"]], exprs87411_Entrez_test, 
                         exprs105777_Entrez_test, exprs126870_Entrez_test)
names(Entrez_exprs_test) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                             "E4_1", "E4_2")

HGNC_exprs_trainval = list(Loaded_exprs_HGNC_trainval[["C1"]], exprs32603_HGNC_trainval, 
                           Loaded_exprs_HGNC_trainval[["C3"]], exprs20181_HGNC_trainval, 
                           exprs55374_HGNC_trainval, exprs59515_HGNC_trainval, 
                  Loaded_exprs_HGNC_trainval[["E2"]], exprs87411_HGNC_trainval, 
                  exprs105777_HGNC_trainval, exprs126870_HGNC_trainval)
names(HGNC_exprs_trainval) = names(Entrez_exprs_trainval)

HGNC_exprs_test = list(Loaded_exprs_HGNC_test[["C1"]], exprs32603_HGNC_test, 
                       Loaded_exprs_HGNC_test[["C3"]], exprs20181_HGNC_test, 
                       exprs55374_HGNC_test, exprs59515_HGNC_test, 
                       Loaded_exprs_HGNC_test[["E2"]], exprs87411_HGNC_test, 
                       exprs105777_HGNC_test, exprs126870_HGNC_test)
names(HGNC_exprs_test) = names(Entrez_exprs_test)

dannot_Entrez_trainval = list(Entrez_dannot_trainval[["C1"]], dannot32603_Entrez_trainval, 
                              Entrez_dannot_trainval[["C3"]], dannot20181_Entrez_trainval, 
                              dannot55374_Entrez_trainval, dannot59515_Entrez_trainval, 
                     Entrez_dannot_trainval[["E2"]], dannot87411_Entrez_trainval, 
                     dannot105777_Entrez_trainval, dannot126870_Entrez_trainval)
names(dannot_Entrez_trainval) = names(Entrez_exprs_trainval)

dannot_Entrez_test = list(Entrez_dannot_test[["C1"]], dannot32603_Entrez_test, 
                          Entrez_dannot_test[["C3"]], dannot20181_Entrez_test, 
                          dannot55374_Entrez_test, dannot59515_Entrez_test, 
                          Entrez_dannot_test[["E2"]], dannot87411_Entrez_test, 
                          dannot105777_Entrez_test, dannot126870_Entrez_test)
names(dannot_Entrez_test) = names(Entrez_exprs_test)

dannot_HGNC_trainval = list(HGNC_dannot_trainval[["C1"]], dannot32603_HGNC_trainval, 
                            HGNC_dannot_trainval[["C3"]], dannot20181_HGNC_trainval, 
                            dannot55374_HGNC_trainval, dannot59515_HGNC_trainval, 
                   HGNC_dannot_trainval[["E2"]], dannot87411_HGNC_trainval, 
                   dannot105777_HGNC_trainval, dannot126870_HGNC_trainval)
names(dannot_HGNC_trainval) = names(Entrez_exprs_trainval)

dannot_HGNC_test = list(HGNC_dannot_test[["C1"]], dannot32603_HGNC_test, 
                        HGNC_dannot_test[["C3"]], dannot20181_HGNC_test, 
                        dannot55374_HGNC_test, dannot59515_HGNC_test, 
                        HGNC_dannot_test[["E2"]], dannot87411_HGNC_test, 
                        dannot105777_HGNC_test, dannot126870_HGNC_test)
names(dannot_HGNC_test) = names(Entrez_exprs_test)

rm(exprs20181_Entrez_trainval, exprs32603_Entrez_trainval, exprs55374_Entrez_trainval,
   exprs59515_Entrez_trainval, exprs87411_Entrez_trainval, Loaded_exprs_Entrez_trainval, 
   exprs20181_HGNC_trainval, exprs32603_HGNC_trainval, exprs55374_HGNC_trainval, 
   exprs59515_HGNC_trainval, exprs87411_HGNC_trainval, Loaded_exprs_HGNC_trainval, 
   Loaded_sets_trainval, Entrez_dannot_trainval, HGNC_dannot_trainval, 
   dannot20181_HGNC_trainval, dannot32603_HGNC_trainval, dannot55374_HGNC_trainval, 
   dannot59515_HGNC_trainval, dannot87411_HGNC_trainval, dannot20181_Entrez_trainval, 
   dannot32603_Entrez_trainval, dannot55374_Entrez_trainval, dannot59515_Entrez_trainval, 
   dannot87411_Entrez_trainval, pam50, exprs105777_Entrez_trainval, exprs105777_HGNC_trainval, 
   exprs126870_Entrez_trainval, exprs126870_HGNC_trainval, dannot105777_Entrez_trainval, 
   dannot126870_Entrez_trainval, dannot105777_HGNC_trainval, dannot126870_HGNC_trainval,
   exprs20181_Entrez_test, exprs32603_Entrez_test, exprs55374_Entrez_test,
   exprs59515_Entrez_test, exprs87411_Entrez_test, Loaded_exprs_Entrez_test, 
   exprs20181_HGNC_test, exprs32603_HGNC_test, exprs55374_HGNC_test, 
   exprs59515_HGNC_test, exprs87411_HGNC_test, Loaded_exprs_HGNC_test, 
   Loaded_sets_test, Entrez_dannot_test, HGNC_dannot_test, 
   dannot20181_HGNC_test, dannot32603_HGNC_test, dannot55374_HGNC_test, 
   dannot59515_HGNC_test, dannot87411_HGNC_test, dannot20181_Entrez_test, 
   dannot32603_Entrez_test, dannot55374_Entrez_test, dannot59515_Entrez_test, 
   dannot87411_Entrez_test, exprs105777_Entrez_test, exprs105777_HGNC_test, 
   exprs126870_Entrez_test, exprs126870_HGNC_test, dannot105777_Entrez_test, 
   dannot126870_Entrez_test, dannot105777_HGNC_test, dannot126870_HGNC_test)
gc()

# Genefu phenotypic annotation #####

Genefu_sets_trainval = list()
for (i in 1:length(Entrez_exprs_trainval)){
  genefu_set = Genefu_predictions_2(expressionDF = Entrez_exprs_trainval[[i]],
                                  HGNC = HGNC_exprs_trainval[[i]],
                                  dannot_Entrez = dannot_Entrez_trainval[[i]],
                                  dannot_HGNC = dannot_HGNC_trainval[[i]])
  Genefu_sets_trainval[[i]] = genefu_set
}

rm(claudinLowData, genefu_set, pam50.robust, scmod1.robust, sig.gene70, i); gc()
names(Genefu_sets_trainval) = names(GEOsets)

Genefu_sets_test = list()
for (i in 1:length(Entrez_exprs_test)){
  genefu_set = Genefu_predictions_2(expressionDF = Entrez_exprs_test[[i]],
                                  HGNC = HGNC_exprs_test[[i]],
                                  dannot_Entrez = dannot_Entrez_test[[i]],
                                  dannot_HGNC = dannot_HGNC_test[[i]])
  Genefu_sets_test[[i]] = genefu_set
}

rm(claudinLowData, genefu_set, pam50.robust, scmod1.robust, sig.gene70, i); gc()
names(Genefu_sets_test) = names(GEOsets)

Genefu_megaset_trainval = rbind(Genefu_sets_trainval[["C1"]], Genefu_sets_trainval[["C2"]], 
                                Genefu_sets_trainval[["C3"]], Genefu_sets_trainval[["E1_1"]], 
                                Genefu_sets_trainval[["E1_2"]], Genefu_sets_trainval[["E1_3"]],
                                Genefu_sets_trainval[["E2"]], Genefu_sets_trainval[["E3"]],
                                Genefu_sets_trainval[["E4_1"]], Genefu_sets_trainval[["E4_2"]])

Genefu_megaset_test = rbind(Genefu_sets_test[["C1"]], Genefu_sets_test[["C2"]], 
                            Genefu_sets_test[["C3"]], Genefu_sets_test[["E1_1"]], 
                            Genefu_sets_test[["E1_2"]], Genefu_sets_test[["E1_3"]],
                            Genefu_sets_test[["E2"]], Genefu_sets_test[["E3"]],
                            Genefu_sets_test[["E4_1"]], Genefu_sets_test[["E4_2"]])

gc()

Pheno_trainval = rbind(train_set, validation_set) %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")
Pheno_train = train_set %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")
Pheno_val = validation_set %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")

# Manually append test set with the dropped samples prior to joining:
rogues = setdiff(Genefu_megaset_test$Sample.ID, test_set$Sample.ID)
one = c("5105", rogues[1], "C2.1001", "C2", "Chemotherapy", "1", "1", "0",
        "RCB, cutoff: >1", "Non_responder", "0", "Pre-treatment", "T1", "GPL14668",
        "Agilent")
two = c("5145", rogues[2], "E1.1.10", "E1_1", "Endocrine_treatment", "1", "0", "1",
        "US_tumour_volume", "Responder", "1", "Pre-treatment", "T1", "GPL96",
        "Affymetrix")
three = c("6145", rogues[3], "E1.2.182", "E1_2", "Endocrine_treatment", "1", "0", "1",
          "US_tumour_volume", "Responder", "1", "Pre-treatment", "T1", "GPL10558",
          "Illumina")
four = c("4145", rogues[4], "E1.3.150", "E1_3", "Endocrine_treatment", "1", "0", "1",
         "US_tumour_volume", "Non_responder", "0", "Pre-treatment", "T1", "GPL10558",
         "Illumina")
five = c("4195", rogues[5], "E3.16723", "E3", "Endocrine_treatment", "1", "0", "1",
         "pCR (+Ki67)", "Non_responder", "0", "Pre-treatment", "T1", "GPL6480",
         "Agilent")
six = c("4159", rogues[6], "E4.1", "E4_1", "Endocrine_treatment", "1", "0", "1",
        "Ki67 60% expr. decr.", "Responder", "1", "Pre-treatment", "T1", "GPL10558",
        "Illumina")
seven = c("4165", rogues[7], "E4.174", "E4_2", "Endocrine_treatment", "1", "0", "1",
          "Ki67 60% expr. decr.", "Non_responder", "0", "Pre-treatment", "T1", "GPL10558",
          "Illumina")

test_set = rbind(test_set, one, two, three, four, five, six, seven)
Pheno = rbind(Pheno, one, two, three, four, five, six, seven)
Pheno_test = test_set %>%
  inner_join(Genefu_megaset_test, by = "Sample.ID")
Pheno_full = rbind(Pheno_trainval, Pheno_test)
Pheno_full = Pheno_full[order(Pheno_full$Dataset),]
rm(Genefu_megaset_trainval, Genefu_megaset_test, Genefu_sets_trainval, Genefu_sets_test,
   dannot_Entrez_trainval, dannot_Entrez_test, dannot_HGNC_trainval,
   dannot_HGNC_test, HGNC_exprs_trainval, HGNC_exprs_test, one, two, three, four, five,
   six, seven, rogues); gc()

Pairwise = inner_join(Pheno[which(Pheno$Timepoint_coded == "T1"),],
                      Pheno[which(Pheno$Timepoint_coded == "T2"),], 
                      by = c("Patient.ID", "Treatment", "Treatment_status",
                             "Chemo_status", "Endo_status", "Response_type", "Response",
                             "Response_coded", "Platform", "Platform_comp")) %>%
  dplyr::select(-Timepoint.y, -Timepoint_coded.y, -Timepoint.x, -Timepoint_coded.x)

# Overview of samples, patients and pre-/on-pairs
rbind(table(Pheno$Dataset), cbind(length(unique(Pheno$Patient.ID[Pheno$Dataset == "C1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "C2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "C3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E4_1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E4_2"]))), 
      table(Pairwise$Dataset.x), table(Pairwise$Dataset.y))

# Taking into account the fact that the last two datasets have some overlap, we have:
# 1610 (out of original 1858) samples, 824 patients and 545 pre-/on-pairs

phenowb = createWorkbook()
addWorksheet(phenowb, "Full")
writeData(phenowb, "Full", Pheno_full)
addWorksheet(phenowb, "Train_val")
writeData(phenowb, "Train_val", Pheno_trainval)
addWorksheet(phenowb, "Train")
writeData(phenowb, "Train", Pheno_train)
addWorksheet(phenowb, "Validation")
writeData(phenowb, "Validation", Pheno_val)
addWorksheet(phenowb, "Test")
writeData(phenowb, "Test", Pheno_test)
saveWorkbook(phenowb, "data/Output sets/Pheno.xlsx", overwrite = TRUE)


detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)
# QC plots #####
# Preparation

library(reshape2)
library(pheatmap)
library(ggplot2)

exprs = inner_join(Entrez_exprs_trainval[[1]], Entrez_exprs_trainval[[2]], 
                   by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[3]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[4]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[5]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[6]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[7]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[8]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[9]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[10]], by = "EntrezGene.ID")

exprs = na.omit(exprs)
Pheno_exprs = Pheno_trainval
check_exprs = as.matrix(exprs[, 2:ncol(exprs)])
check_exprs = check_exprs[, Pheno_exprs$Sample.ID]
rownames(check_exprs) = exprs$EntrezGene.ID

plot_matrix = as.matrix(exprs[,2:ncol(exprs)])
plot_matrix = plot_matrix[, Pheno_exprs$Sample.ID]
rownames(plot_matrix) = exprs$EntrezGene.ID
mds = plotMDS(plot_matrix)
pca = data.frame(cbind(mds$x, mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                       Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                       as.character(Pheno_exprs$Timepoint_coded),
                       Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                       as.character(Pheno_exprs$IC10), 
                       as.character(Pheno_exprs$Mammaprint_risk), 
                       as.character(Pheno_exprs$rorS_risk)))
pca$Mammaprint = pca$X11
pca$Mammaprint[which(pca$X11=="1")] = "Risk"
pca$Mammaprint[which(pca$X11!="1")] = "No risk"
pca = pca %>% dplyr::select(-X11)
colnames(pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                  "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  pca[,i] = factor(pca[,i], levels = unique(pca[,i]), labels = unique(pca[,i]))
}

class(pca$X1) = "numeric"
class(pca$X2) = "numeric"

# Multidimensional scaling plots #####
gpca = ggplot(pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS22 (", round(100*mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
gpca
ggsave(filename = "MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# RBE
RBE_design = model.matrix(~0 + Pheno_exprs$Response + 
                            Pheno_exprs$Timepoint_coded +
                            Pheno_exprs$pam50)[,c(1:3, 6:9)]
colnames(RBE_design) = c("Non_responder", "Responder", "T2",
                         "HER2", "LumB", "LumA", "Normal")
rownames(RBE_design) = Pheno_exprs$Sample.ID
cm_RBE = makeContrasts(RespvsNonresp = Responder - Non_responder,
                       levels = RBE_design)
seq_batch = Pheno_exprs$Dataset
RBE_exprs = removeBatchEffect(check_exprs, batch = seq_batch, design = RBE_design)

RBE_mds = plotMDS(RBE_exprs)
RBE_pca = data.frame(cbind(RBE_mds$x, RBE_mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                           Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                           as.character(Pheno_exprs$Timepoint_coded),
                           Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                           as.character(Pheno_exprs$IC10), 
                           as.character(Pheno_exprs$Mammaprint_risk), 
                           as.character(Pheno_exprs$rorS_risk)))
RBE_pca$Mammaprint = RBE_pca$X11
RBE_pca$Mammaprint[which(RBE_pca$X11=="1")] = "Risk"
RBE_pca$Mammaprint[which(RBE_pca$X11!="1")] = "No risk"
RBE_pca = RBE_pca %>% dplyr::select(-X11)
colnames(RBE_pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                      "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  RBE_pca[,i] = factor(RBE_pca[,i], levels = unique(RBE_pca[,i]), labels = unique(RBE_pca[,i]))
}

class(RBE_pca$X1) = "numeric"
class(RBE_pca$X2) = "numeric"

rbe_gpca = ggplot(RBE_pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "RBE - Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*RBE_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS2 (", round(100*RBE_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
rbe_gpca
ggsave(filename = "RBE_MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# z-score normalisation ( https:://www.biostars.org/p/283083/ )
z = list()
means = list()
sd = list()

for(i in 1:length(Entrez_exprs_trainval)){
  df = as.data.frame(Entrez_exprs_trainval[[i]][,2:ncol(Entrez_exprs_trainval[[i]])])
  t = as.data.frame(t(df))
  means[[i]] = sapply(t, function(t) mean(t, na.rm = T))
  sd[[i]] = sapply(t, function(t) sd(t, na.rm = T))
  z_t = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
  z[[i]] = as.matrix(t(z_t))
  rownames(z[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
  colnames(z[[i]]) = colnames(df)
  z[[i]] = as.data.frame(z[[i]])
  z[[i]]$EntrezGene.ID = Entrez_exprs_trainval[[i]]$EntrezGene.ID
  rm(t, z_t, df)
}

for (i in 1:length(means)) names(means[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
for (i in 1:length(sd)) names(sd[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
names(z) = names(means) = names(sd) = names(Entrez_exprs_trainval)

# To apply normalisation consistently on the test set we will use the saved 
# means and standard deviations in the means and sd lists that we created:

z_test = list()
for (i in 1:length(Entrez_exprs_test)){
  df = as.data.frame(Entrez_exprs_test[[i]][,2:ncol(Entrez_exprs_test[[i]])])
  t = as.data.frame(t(df))
  z_t = t
  for (j in 1:nrow(t)){
    for (k in 1:ncol(t)){
      z_t[j, k] = (z_t[j, k] - means[[i]][k])/sd[[i]][k]
    }
  }
  z_test[[i]] = as.matrix(t(z_t))
  rownames(z_test[[i]]) = Entrez_exprs_test[[i]]$EntrezGene.ID
  colnames(z_test[[i]]) = colnames(df)
  z_test[[i]] = as.data.frame(z_test[[i]])
  z_test[[i]]$EntrezGene.ID = Entrez_exprs_test[[i]]$EntrezGene.ID
  rm(t, z_t, df)
}

names(z_test) = names(z)

# Joining in one expression matrix
# Test
z_test_exprs = z_test[[1]] %>% inner_join(z_test[[2]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[3]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[4]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[5]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[6]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[7]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[8]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[9]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[10]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

z_test_exprs = na.omit(z_test_exprs)
ids = z_test_exprs$EntrezGene.ID
z_test_exprs = as.matrix(z_test_exprs %>% dplyr::select(-EntrezGene.ID))
z_test_exprs = z_test_exprs[, Pheno_test$Sample.ID]
rownames(z_test_exprs) = ids
z_test_exprs = z_test_exprs[rowSums(is.na(z_test_exprs)) != ncol(z_test_exprs), ]

# Training and Validation
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(z[[4]], by = "EntrezGene.ID") %>%
  inner_join(z[[5]], by = "EntrezGene.ID") %>%
  inner_join(z[[6]], by = "EntrezGene.ID") %>%
  inner_join(z[[7]], by = "EntrezGene.ID") %>%
  inner_join(z[[8]], by = "EntrezGene.ID") %>%
  inner_join(z[[9]], by = "EntrezGene.ID") %>%
  inner_join(z[[10]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

z_exprs = na.omit(z_exprs)
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID))
z_exprs = z_exprs[, Pheno_exprs$Sample.ID]
rownames(z_exprs) = rownames(check_exprs)
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
KBZ_mds = plotMDS(z_exprs)
KBZ_pca = data.frame(cbind(KBZ_mds$x, KBZ_mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                           Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                           as.character(Pheno_exprs$Timepoint_coded),
                           Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                           as.character(Pheno_exprs$IC10), 
                           as.character(Pheno_exprs$Mammaprint_risk), 
                           as.character(Pheno_exprs$rorS_risk)))
KBZ_pca$Mammaprint = KBZ_pca$X11
KBZ_pca$Mammaprint[which(KBZ_pca$X11=="1")] = "Risk"
KBZ_pca$Mammaprint[which(KBZ_pca$X11!="1")] = "No risk"
KBZ_pca = KBZ_pca %>% dplyr::select(-X11)
colnames(KBZ_pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                      "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  KBZ_pca[,i] = factor(KBZ_pca[,i], levels = unique(KBZ_pca[,i]), labels = unique(KBZ_pca[,i]))
}

class(KBZ_pca$X1) = "numeric"
class(KBZ_pca$X2) = "numeric"

z_pca = ggplot(KBZ_pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Normalised samples - Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*KBZ_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS2 (", round(100*KBZ_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
z_pca
ggsave(filename = "z_MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Expression boxplots - Not helpful due to large number of samples #####
eset = exprs[,2:ncol(exprs)]
png("QC/Boxplot.png", width = 1920, height = 1080)
ggplot(melt(eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C1",])), 
                        rep("orange", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C2",])),
                        rep("red4", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C3",])),
                        rep("darkseagreen", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_1",])),
                        rep("beige", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_2",])),
                        rep("blue4", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_3",])),
                        rep("deeppink", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E2",])),
                        rep("purple", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E3",])),
                        rep("grey", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E4_1",])),
                        rep("darkgreen", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E4_2",]))), 
               outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression",
       x = "\nSamples",
       y = "Expression\n")
dev.off()

# Heatmaps #####
save_pheatmap_png <- function(x, filename, width = 174, height = 120.8333, res = 650,
                              units = "mm") {
  png(filename, width = width, height = height, res = res, units = units, 
      compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_tiff <- function(x, filename, width = 174, height = 120.8333, res = 650,
                               units = "mm") {
  tiff(filename, width = width, height = height, res = res, units = units, 
       compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


annotation_for_heatmap = pca[, c(3, 5:12)]
rownames(annotation_for_heatmap) = pca$Sample.ID

dists = as.matrix(dist(t(plot_matrix), method = "manhattan"))

rownames(dists) = pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Response = c(Responder = "dodgerblue4", Non_responder = "deeppink4"),
  Dataset = c(C1 = "cyan", C2 = "orange", C3 = "red4", E1_1 = "darkseagreen",
              E1_2 = "beige", E1_3 = "blue4", E2 = "deeppink",
              E3 = "purple", E4_1 = "grey", E4_2 = "darkgreen"),
  Timepoint = c(T1 = "goldenrod2", T2 = "purple4"),
  Treatment = c(Chemotherapy = "paleturquoise3", Endocrine_treatment = "sienna2"),
  scmod1 = c(`ER-/HER2-` = "red4", `ER+/HER2- High Prolif` = "skyblue",
             `ER+/HER2- Low Prolif` = "darkblue", `HER2+` = "violet"),
  pam50 = c(Basal = "red4", LumB = "skyblue",
            LumA = "darkblue", Her2 = "purple4", Normal = "lightgreen"),
  IC10 = c(iC1 = "pink3", iC2 = "red", iC3 = "red4", iC4 = "beige", iC5 = "aliceblue",
           iC6 = "cadetblue3", iC7 = "darkmagenta", iC8 = "hotpink4", iC9 = "lightpink2",
           iC10 = "mistyrose1"),
  rorS = c(High = "red4", Intermediate = "orange", Low = "chartreuse4"),
  Mammaprint = c(Risk = "red3", `No risk` = "green")
)

heatmap = pheatmap(t(dists), col = (hmcol), 
                   annotation_col = annotation_for_heatmap,
                   annotation_colors = ann_colors,
                   cluster_cols = T,
                   treeheight_col =  0,
                   legend = TRUE,
                   show_rownames = F,
                   show_colnames = F,
                   fontsize = 3.5,
                   legend_breaks = c(min(dists, na.rm = TRUE), 
                                     max(dists, na.rm = TRUE)), 
                   legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(heatmap, "QC/Original_heatmap.tiff")

# RBE
RBE_annotation_for_heatmap = RBE_pca[, c(3, 5:12)]
rownames(RBE_annotation_for_heatmap) = RBE_pca$Sample.ID

rbe_dists = as.matrix(dist(t(RBE_exprs), method = "manhattan"))

rownames(rbe_dists) = RBE_pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(rbe_dists) <- NULL
diag(rbe_dists) <- NA

RBE_heatmap = pheatmap(t(rbe_dists), col = (hmcol), 
                       annotation_col = RBE_annotation_for_heatmap,
                       annotation_colors = ann_colors,
                       cluster_cols = T,
                       treeheight_col =  0,
                       legend = TRUE,
                       show_colnames = F,
                       show_rownames = F,
                       fontsize = 3.5,
                       legend_breaks = c(min(rbe_dists, na.rm = TRUE), 
                                         max(rbe_dists, na.rm = TRUE)), 
                       legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(RBE_heatmap, "QC/RBE_heatmap.tiff")

# KBZ
KBZ_annotation_for_heatmap = KBZ_pca[, c(3, 5:12)]
rownames(KBZ_annotation_for_heatmap) = KBZ_pca$Sample.ID

z_dists = as.matrix(dist(t(z_exprs), method = "manhattan"))

rownames(z_dists) = KBZ_pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

KBZ_heatmap = pheatmap(t(z_dists), col = (hmcol), 
                       annotation_col = KBZ_annotation_for_heatmap,
                       annotation_colors = ann_colors,
                       cluster_cols = T,
                       treeheight_col =  0,
                       legend = TRUE,
                       show_colnames = F,
                       show_rownames = F,
                       fontsize = 3.5,
                       legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                         max(z_dists, na.rm = TRUE)), 
                       legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(KBZ_heatmap, "QC/z_heatmap.tiff")

# 3D plots
library(rgl)
library(plot3D)
library(plot3Drgl)

exprs_for_join = exprs
str_sub(exprs_for_join$EntrezGene.ID, 0, 0) = "G_"
rownames(exprs_for_join) = exprs_for_join$EntrezGene.ID
exprs_for_join = as.data.frame(t(exprs_for_join))
exprs_for_join = exprs_for_join[2:nrow(exprs_for_join),]
exprs_for_join$Sample.ID = rownames(exprs_for_join)
masternormal = pca %>% inner_join(exprs_for_join, by =  "Sample.ID")

normaldist <- dist(masternormal[,13:ncol(masternormal)])
normalmds <- cmdscale(normaldist, k = 3)
PC1 <- normalmds[,1]
PC2 <- normalmds[,2]
PC3 <- normalmds[,3]

#Plot
scatter3D(PC1, PC2, PC3, phi = 10, bty = "u", theta = 45, width = 1920, height = 1080,
          main = "3D PCA - raw data", xlab = "PC1", ylab = "PC2", zlab = "PC3",
          ticktype = "detailed", pch = 20, alpha = 0.7,
          colvar = as.integer(masternormal$Dataset),
          col = c("cyan", "orange", "red4", "darkseagreen",
                  "beige", "blue4", "deeppink",
                  "purple", "grey", "darkgreen"), 
          colkey = list(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                        addlines = TRUE, length = 1, width = 1,
                        labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3",
                                   "E2", "E3", "E4_1", "E4_2")), clab = "Dataset", cex = 1)
plotrgl()

# KBZ #
z_exprs_for_join = as.data.frame(t(z_exprs))
str_sub(colnames(z_exprs_for_join), 0, 0) = "G_"
z_exprs_for_join$Sample.ID = rownames(z_exprs_for_join)
KBZ_master = KBZ_pca %>% inner_join(z_exprs_for_join, by =  "Sample.ID")

KBZdist <- dist(KBZ_master[,13:ncol(KBZ_master)])
KBZmds <- cmdscale(KBZdist, k = 3)
KBZPC1 <- KBZmds[,1]
KBZPC2 <- KBZmds[,2]
KBZPC3 <- KBZmds[,3]

# Plot
scatter3D(KBZPC1, KBZPC2, KBZPC3, phi = 10, bty = "u", theta = 45, width = 1920, height = 1080,
          main = "3D PCA - z-scores", xlab = "PC1", ylab = "PC2", zlab = "PC3",
          ticktype = "detailed", pch = 20, alpha = 0.7,
          colvar = as.integer(KBZ_master$Dataset),
          col = c("cyan", "orange", "red4", "darkseagreen",
                  "beige", "blue4", "deeppink",
                  "purple", "grey", "darkgreen"), 
          colkey = list(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                        addlines = TRUE, length = 1, width = 1,
                        labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3",
                                   "E2", "E3", "E4_1", "E4_2")), clab = "Dataset", cex = 1)
plotrgl()

rm(i, j, k)

##### Plotly #####
library(plotly)
library(data.table)
library(colorspace)

Pheno_sunburst = Pheno_full %>% dplyr::select(Treatment, Dataset, Response, pam50)

Pheno_sunburst = Pheno_sunburst %>%
  group_by(Treatment, Dataset, Response, pam50) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

as.sunburstDF <- function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF <- names(DF)
  
  if(is.data.table(DF)){
    DT <- copy(DF)
  } else {
    DT <- data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT <- names(DT)
  hierarchy_columns <- setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList <- list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns <- colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT <- unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT <- DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <- rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns <- setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}

coloring  = data.frame(stringsAsFactors = FALSE,
                       colors = tolower(gplots::col2hex(c("paleturquoise3", "sienna2", "paleturquoise1",
                                                          "paleturquoise1", "paleturquoise1", "sienna1", 
                                                          "sienna1", "sienna1", "sienna1",
                                                          "sienna1", "sienna1", "sienna1", "deeppink4", "dodgerblue4", 
                                                          "red4", "violet", "darkblue", "skyblue", "lightgreen"))),
                       labels = c("Chemotherapy", "Endocrine_treatment", "C1", "C2", "C3", 
                                  "E1_1", "E1_2", "E1_3", "E2", "E3", "E4_1", "E4_2",
                                  "Non_responder", "Responder", "Basal", "Her2", "LumA", 
                                  "LumB", "Normal"))

sunburstDF <- as.sunburstDF(Pheno_sunburst, value_column = "Counts", add_root = FALSE) %>%
  inner_join(coloring, by = "labels")

pie = plot_ly() %>%
  add_trace(ids = sunburstDF$ids, labels= sunburstDF$labels, parents = sunburstDF$parents, 
            values= sunburstDF$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 4,
            marker = list(colors = sunburstDF$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie

# Session Info #####
sessionInfo()
=======
# Libraries #####
library(limma)
library(org.Hs.eg.db)
library(dplyr)
library(openxlsx)
library(GEOquery)
library(genefu)
library(stringr)
library(tidyr)
library(EnhancedVolcano)
library(readr)
# library(knitr)
# library(rJava)
readr::local_edition(1) # required for v>=2.60.0 of GEOquery

# Loading data #####

# Downloading the GEO objects for all studies of interest
tvt_series = c("GSE32603", "GSE20181", "GSE55374", "GSE59515", "GSE87411",
               "GSE105777", "GSE126870")

ext_val_series = c("GSE18728", "GSE119262", "GSE111563")

GEOsets_pre = list()
for (i in 1:length(tvt_series)){
  GEOsets_pre[[i]] = unlist(getGEO(tvt_series[i])[[1]])
}
names(GEOsets_pre) = tvt_series; rm(tvt_series, i); gc()

Ext_val = list()
for (i in 1:length(ext_val_series)){
  Ext_val[[i]] = unlist(getGEO(ext_val_series[i])[[1]])
}
names(Ext_val) = ext_val_series; rm(ext_val_series, i); gc()

# Expression matrices for the Bownes, Dunbier and Park studies will be loaded 
# as local files. We stored the downloaded objects in "data/Expression Matrices"
# Loading the remaining Expression Matrices
detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)
Bownes = read.csv("data/Expression Matrices/matrixC1.csv") # available in the andysims lab
Dunbier = read.csv("data/Expression Matrices/matrix_E2.csv") # not available online - provided by Anita Dunbier to andysims lab
Park = read.csv("data/Expression Matrices/matrixC3.csv") # GEOquery fails in retrieving the data
# Removing some unidentified probes
Park = Park %>% dplyr::slice(1:16692,16706:30373)

GEOsets = list(Bownes, GEOsets_pre[[1]], Park, GEOsets_pre[[2]], GEOsets_pre[[3]], 
               GEOsets_pre[[4]], Dunbier, GEOsets_pre[[5]], GEOsets_pre[[6]],
               GEOsets_pre[[7]])
names(GEOsets) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                   "E4_1", "E4_2")

rm(GEOsets_pre, Bownes, Park, Dunbier)

# Phenotypic annotation #####
# Entrez_pheno (HGNC_pheno is the same essentially)
Entrez_pheno = GEOsets
for (i in 1:length(GEOsets)){
  Entrez_pheno[[i]] = as.data.frame(colnames(GEOsets[[i]][,2:ncol(GEOsets[[i]])]))
  colnames(Entrez_pheno[[i]]) = "Sample.ID"
  Entrez_pheno[[i]]$Dataset = names(Entrez_pheno[i])
  if(i <= 3){Entrez_pheno[[i]]$Treatment = "Chemotherapy"} else 
  {Entrez_pheno[[i]]$Treatment = "Endocrine_treatment"}
}

# Starting from C1 - info taken from publication and exprs sample names
Entrez_pheno[[1]]$Timepoint_coded = str_sub(Entrez_pheno[[1]]$Sample.ID, -2, -1)
Entrez_pheno[[1]]$Timepoint = as.character(Entrez_pheno[[1]]$Timepoint_coded)
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T2"] = "2 weeks"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T3"] = "1.5 months"
Entrez_pheno[[1]]$Timepoint[Entrez_pheno[[1]]$Timepoint_coded == "T4"] = "Surgery"
Entrez_pheno[[1]]$Timepoint_coded[Entrez_pheno[[1]]$Timepoint == "1.5 months"] = "T2.5"
Entrez_pheno[[1]]$Timepoint_coded[Entrez_pheno[[1]]$Timepoint == "Surgery"] = "T3"
Entrez_pheno[[1]]$Patient.ID = str_sub(Entrez_pheno[[1]]$Sample.ID, 1,
                                       nchar(Entrez_pheno[[1]]$Sample.ID) - 3)
str_sub(Entrez_pheno[[1]]$Patient.ID, 1, 3) <- "C1."
Entrez_pheno[[1]]$Response_type = "pCR"
Entrez_pheno[[1]]$Response = NA

# C2
# times
pdata2 = pData(GEOsets[["C2"]])
times = data.frame(rownames(pdata2),pdata2$characteristics_ch1.3)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times$Timepoint_coded = str_sub(times$Timepoint_coded, -2, -1)
times[which(times$Timepoint_coded == "TS"),2] = "T3"
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% inner_join(times, by = "Sample.ID")
Entrez_pheno[[2]]$Timepoint = as.character(Entrez_pheno[[2]]$Timepoint_coded)
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T2"] = "24-96 hours"
Entrez_pheno[[2]]$Timepoint[Entrez_pheno[[2]]$Timepoint_coded == "T3"] = "Surgery"

# patient ID's
patient_ids = data.frame(rownames(pdata2), pdata2$characteristics_ch1)
patient_ids = patient_ids %>%
  dplyr::rename(Sample.ID=colnames(patient_ids)[1], 
                Patient.ID=colnames(patient_ids)[2])
patient_ids$Patient.ID = str_sub(patient_ids$Patient.ID, -4, -1)
str_sub(patient_ids$Patient.ID, 0, 0) = "C2."
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% left_join(patient_ids, by = "Sample.ID")

# Response
Entrez_pheno[[2]]$Response_type = "RCB, cutoff: >1"
response = data.frame(rownames(pdata2), pdata2$characteristics_ch1.9)
response = response %>%
  dplyr::rename(Sample.ID=colnames(response)[1], 
                Response=colnames(response)[2])
response$Response = str_sub(response$Response, 11, -1)
response$Response = as.numeric(response$Response)
response[which(response$Response == 3),2] = "Non_responder"
response[which(response$Response == 2),2] = "Non_responder"
response[which(response$Response == 1),2] = "Responder"
response[which(response$Response == 0),2] = "Responder"
Entrez_pheno[[2]] = Entrez_pheno[[2]] %>% inner_join(response, by = "Sample.ID")

rm(pdata2, times, patient_ids, response)

# C3
# times
library(GEOquery)
readr::local_edition(1)
pdata3 = pData(getGEO('GSE123845',GSEMatrix=TRUE)[[1]])
times = data.frame(pdata3$title, pdata3$characteristics_ch1.54)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times$Timepoint_coded = str_sub(times$Timepoint_coded, -2, -1)
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(times, by = "Sample.ID")
Entrez_pheno[[3]]$Timepoint = as.character(Entrez_pheno[[3]]$Timepoint_coded)
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T1"] = "Pre-treatment"
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T2"] = "3 weeks"
Entrez_pheno[[3]]$Timepoint[Entrez_pheno[[3]]$Timepoint_coded == "T3"] = "6 months/Surgery"

# patient ID's
patient_ids = data.frame(pdata3$title)
patient_ids$Patient.ID = as.character(patient_ids$pdata3.title)
patient_ids = patient_ids %>%
  dplyr::rename(Sample.ID=colnames(patient_ids)[1])
patient_ids$Patient.ID = str_sub(patient_ids$Patient.ID, 4, -3)
str_sub(patient_ids$Patient.ID, 0, 0) = "C3."
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(patient_ids, by = "Sample.ID")

# Response
response = data.frame(pdata3$title, pdata3$characteristics_ch1.35)
response = response %>%
  dplyr::rename(Sample.ID=colnames(response)[1], 
                Response=colnames(response)[2])
response$Response_type = "pCR/RD"
response[which(response$Response == "pcr_status: 0"),2] = "Responder"
response[which(response$Response == "pcr_status: 1"),2] = "Non_responder"
Entrez_pheno[[3]] = Entrez_pheno[[3]] %>% inner_join(response, by = "Sample.ID")

rm(pdata3, times, patient_ids, response)
detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)

# E1_1
# times & patient ID's & response
pdata4 = pData(GEOsets[["E1_1"]])
times = data.frame(rownames(pdata4), pdata4$title)
times = times %>%
  dplyr::rename(Sample.ID=colnames(times)[1], Timepoint_coded=colnames(times)[2])
times = separate(times, Timepoint_coded, into = c("Timepoint_coded",
                                                  "col2", "col3", "col4",
                                                  "Response"),
                 sep = ";")
times$Response[117:176] = times$col2[117:176]
times = times %>% dplyr::select(-col2, -col3, -col4)
times = separate(times, Timepoint_coded, into = c("Timepoint_coded",
                                                  "col2", "col3", "col4"),
                 sep = ":")
times$Response_type = "US_tumour_volume"
times[which(times$Response == " responder"), 6] = "Responder"
times[which(times$Response == " nonresponder"), 6] = "Non_responder"
times[which(times$Response == " not assessable"), 6] = NA

times$Patient.ID = str_sub(times$Timepoint_coded,0,-2)
str_sub(times$Patient.ID, 0, 0) = "E1.1."
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "A"),2] = "T1"
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "B"),2] = "T2"
times[which(substr(times$Timepoint_coded, nchar(times$Timepoint_coded),
                   nchar(times$Timepoint_coded)) == "C"),2] = "T3"
times$Timepoint = as.character(times$Timepoint_coded)
times$Timepoint[times$Timepoint_coded == "T1"] = "Pre-treatment"
times$Timepoint[times$Timepoint_coded == "T2"] = "2 weeks"
times$Timepoint[times$Timepoint_coded == "T3"] = "Surgery"
times = times %>% dplyr::select(Sample.ID, Timepoint_coded, Timepoint,
                                Patient.ID, Response, Response_type)
Entrez_pheno[[4]] = Entrez_pheno[[4]] %>% inner_join(times, by = "Sample.ID")

rm(pdata4, times)

# E1_2
pdata5 = pData(GEOsets[["E1_2"]])
add = pdata5 %>% dplyr::select(geo_accession, `clinical response:ch1`,
                               `subject:ch1`, `treatment:ch1`) %>%
  dplyr::rename(Sample.ID=geo_accession, Response=`clinical response:ch1`,
                Patient.ID=`subject:ch1`, Timepoint=`treatment:ch1`)
add$Response_type = "US_tumour_volume"
str_sub(add$Patient.ID, 0, -4) = "E1.2."
add[which(add$Timepoint == "pre-treatment"), 4] = "Pre-treatment"
add[which(add$Timepoint == "2 wk endocrine therapy"), 4] = "2 weeks"
add[which(add$Timepoint == "3 mo endocrine therapy"), 4] = "3 months/Surgery"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[which(add$Timepoint == "3 months/Surgery")] = "T3"
add$Timepoint_coded[which(add$Timepoint == "2 weeks")] = "T2"
add$Timepoint_coded[which(add$Timepoint == "Pre-treatment")] = "T1"
Entrez_pheno[[5]] = Entrez_pheno[[5]] %>% inner_join(add, by = "Sample.ID")

rm(add, pdata5)

# E1_3
pdata6 = pData(GEOsets[["E1_3"]])
add = pdata6 %>% dplyr::select(geo_accession, `clinical response:ch1`,
                               `patient id:ch1`, `time point:ch1`) %>%
  dplyr::rename(Sample.ID=geo_accession, Response=`clinical response:ch1`,
                Patient.ID=`patient id:ch1`, Timepoint=`time point:ch1`)
add$Response_type = "US_tumour_volume"
str_sub(add$Patient.ID, 0, 0) = "E1.3."
add[which(add$Timepoint == "pre-treatment"), 4] = "Pre-treatment"
add[which(add$Timepoint == "2 wks of letrozole treatment"), 4] = "2 weeks"
add[which(add$Timepoint == "3 months of letrozole treatment"), 4] = "3 months/Surgery"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[which(add$Timepoint == "3 months/Surgery")] = "T3"
add$Timepoint_coded[which(add$Timepoint == "2 weeks")] = "T2"
add$Timepoint_coded[which(add$Timepoint == "Pre-treatment")] = "T1"
Entrez_pheno[[6]] = Entrez_pheno[[6]] %>% inner_join(add, by = "Sample.ID")
Entrez_pheno[[6]]$Response = gsub("Non-responder", "Non_responder", Entrez_pheno[[6]]$Response)

rm(add, pdata6)

# E2
# E2_supp is created in MS Excel from data again available at the Sims lab space
# provided by Anita Dunbier from the Royal Marsden study ("E2")
E2_supp = openxlsx::read.xlsx("data/Supplementary/E2_extras.xlsx", sheet = 1, colNames = T)
E2_supp = E2_supp %>%
  dplyr::select(Sample.ID, Patient.ID, Timepoint,
                `FactorValue[ClinicalInformation-Response]`) %>%
  dplyr::rename(Response = `FactorValue[ClinicalInformation-Response]`)
str_sub(E2_supp$Patient.ID, 0, 0) = "E2."
E2_supp$Timepoint_coded = as.character(E2_supp$Timepoint)
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "16week post-treatment")] = "T3"
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "2 week post-treatment")] = "T2"
E2_supp$Timepoint_coded[which(E2_supp$Timepoint == "pre-treatment")] = "T1"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T1")] = "Pre-treatment"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T2")] = "2 weeks"
E2_supp$Timepoint[which(E2_supp$Timepoint_coded == "T3")] = "16 weeks"
E2_supp$Response = gsub("responder", "Responder", E2_supp$Response)
E2_supp$Response = gsub("non-Responder", "Non_responder", E2_supp$Response)
E2_supp$Response_type = "Ki67 expression"
Entrez_pheno[[7]] = Entrez_pheno[[7]] %>% left_join(E2_supp, by = "Sample.ID")

rm(E2_supp)

# E3
# This set contains reanalysed samples from studies not included in our list
# of studies, so reanalysed samples are kept in this case

pdata8 = pData(GEOsets[["E3"]])
add = pdata8 %>% 
  dplyr::select(`geo_accession`, `source_name_ch2`,`studynum:ch2`, 
                `endocrine therapy response group:ch2`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Timepoint = `source_name_ch2`,
                Patient.ID = `studynum:ch2`, 
                Response = `endocrine therapy response group:ch2`)
str_sub(add$Patient.ID, 0, 0) = "E3."
add[which(add$Response == "resistant"), 4] = "Non_responder"
add[which(add$Response == "sensitive"), 4] = "Responder"
add$Response_type = "pCR (+Ki67)"
add = separate(add, Timepoint, into = c("col1", "col2", "col3",
                                        "Timepoint", "col5"),
               sep = "_")
add = add %>% dplyr::select(-col1, -col2, -col3, -col5)
add[which(add$Timepoint == "BL"), 2] = "Pre-treatment"
add[which(add$Timepoint == "BL Redo"), 2] = "Pre-treatment"
add[which(add$Timepoint == "W2"), 2] = "2 weeks"
add[which(add$Timepoint == "w2"), 2] = "2 weeks"
add[which(add$Timepoint == "w2 Redo"), 2] = "2 weeks"
add[which(add$Timepoint == "M"), 2] = "4 weeks/Mid-treatment"
add$Timepoint_coded = as.character(add$Timepoint)
add$Timepoint_coded[add$Timepoint == "Pre-treatment"] = "T1"
add$Timepoint_coded[add$Timepoint == "2 weeks"] = "T2"
add$Timepoint_coded[add$Timepoint == "4 weeks/Mid-treatment"] = "T2.5"
Entrez_pheno[[8]] = Entrez_pheno[[8]] %>% inner_join(add, by = "Sample.ID")

rm(add, pdata8)

# E4_1
# The set includes additional samples from the POETIC study and the NIT-2 study. 
# Here we will remove the samples from the NIT-2 study that have been reanalysed,
# because these are no-intervention treatment samples from GSE73237. The samples that
# have later been reanalysed in GSE126870 will also be removed.

# Since response is not reported, all of the samples that we are keeping here
# are gonna have missing values on response.

pdata9 = pData(GEOsets[["E4_1"]])
add = pdata9 %>%
  dplyr::select(`geo_accession`, relation, description, `group:ch1`, `timepoint:ch1`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Relation = relation, Desc = description, 
                Treated = `group:ch1`, Timepoint_pre = `timepoint:ch1`) %>%
  dplyr::filter(!grepl("Reanalyzed by: ", Relation)) %>%
  dplyr::filter(!grepl("Reanalysis of: ", Relation)) %>%
  dplyr::select(-Relation, -Treated) # Treated is removed because all subjects were treated
add$Timepoint = add$Timepoint_pre
add$Timepoint[add$Timepoint_pre == "baseline (B)"] = "Pre-treatment"
add$Timepoint[add$Timepoint_pre == "surgery (S)"] = "2 weeks (surgery)"
add$Timepoint_coded = add$Timepoint
add$Timepoint_coded[add$Timepoint == "Pre-treatment"] = "T1"
add$Timepoint_coded[add$Timepoint == "2 weeks (surgery)"] = "T2"
add$Desc = gsub("Treated.", "", add$Desc)
add$Desc = gsub("B", "", add$Desc)
add$Desc = gsub("S", "", add$Desc)
str_sub(add$Desc, 0, 0) = "E4."
add = add %>% dplyr::select(-Timepoint_pre) %>%
  dplyr::rename(Patient.ID = Desc)

# E4_2
# We are going to use the official supplementary file too here, because unlike
# the pheno info deposited in GEO, it contains info on response.
# Only baseline samples (reanalysed or not) are included here. These will be
# matched to the surgery samples from E4_1. However, the info about which study
# each sample came from will be retained to be included as a batch effect

pdata10 = pData(GEOsets[["E4_2"]])
pdata10_ext = openxlsx::read.xlsx("data/Supplementary/E4_2_supp.xlsx", sheet = 5,
                                  startRow = 3)
add1 = pdata10 %>%
  dplyr::select(`geo_accession`, title, `timepoint:ch1`) %>%
  dplyr::rename(Sample.ID = `geo_accession`, Patient.ID = title, Timepoint = `timepoint:ch1`)
add1$Patient.ID = gsub("AI.Treated.", "", add1$Patient.ID)
add1$Patient.ID = gsub(" \\[reanalysis\\]", "", add1$Patient.ID)
add1$Patient.ID = gsub("B", "", add1$Patient.ID)
str_sub(add1$Patient.ID, 0, 0) = "E4."
add1$Timepoint = gsub("baseline \\(B\\)", "Pre-treatment", add1$Timepoint)
add1$Timepoint_coded = "T1"

pdata10_ext = pdata10_ext[, c(2,11)]
colnames(pdata10_ext) = c("Patient.ID", "Response")
add2 = pdata10_ext %>% # Patient and Response columns
  dplyr::filter(!grepl("Control.", Patient.ID))
add2$Response = gsub("responder", "Responder", add2$Response)
add2$Response = gsub("non-Responder", "Non_responder", add2$Response)
add2$Response_type = "Ki67 60% expr. decr."
add2$Patient.ID = gsub("Treated", "E4", add2$Patient.ID)

add3 = add1 %>% left_join(add2, by = "Patient.ID") # for E4_2
add4 = add %>% left_join(add2, by = "Patient.ID") # for E4_1

Entrez_pheno[[9]] = Entrez_pheno[[9]] %>% inner_join(add4, by = "Sample.ID")
Entrez_pheno[[10]] = Entrez_pheno[[10]] %>% inner_join(add3, by = "Sample.ID")

rm(add, add1, add2, add3, add4, pdata10, pdata10_ext, pdata9)

# Adding a coded response column in each dataset
# C1 - using a supplementary data file
C1_Response = openxlsx::read.xlsx("data/Supplementary/C1.Response.xlsx", sheet = 2, colNames = T)
Entrez_pheno[[1]] = Entrez_pheno[[1]] %>% inner_join(C1_Response, by = "Patient.ID")
Entrez_pheno[[1]]$Response[Entrez_pheno[[1]]$Response_coded == 1] = "Responder"
Entrez_pheno[[1]]$Response[Entrez_pheno[[1]]$Response_coded == 0] = "Non_responder"

rm(C1_Response)

# C2
Entrez_pheno[[2]]$Response_coded = NA
Entrez_pheno[[2]]$Response_coded[Entrez_pheno[[2]]$Response == "Responder"] = 1
Entrez_pheno[[2]]$Response_coded[Entrez_pheno[[2]]$Response == "Non_responder"] = 0

# C3
Entrez_pheno[[3]]$Response_coded = NA
Entrez_pheno[[3]]$Response_coded[Entrez_pheno[[3]]$Response == "Responder"] = 1
Entrez_pheno[[3]]$Response_coded[Entrez_pheno[[3]]$Response == "Non_responder"] = 0

# E1_1
Entrez_pheno[[4]]$Response_coded = NA
Entrez_pheno[[4]]$Response_coded[Entrez_pheno[[4]]$Response == "Responder"] = 1
Entrez_pheno[[4]]$Response_coded[Entrez_pheno[[4]]$Response == "Non_responder"] = 0

# E1_2 (all are responders)
Entrez_pheno[[5]]$Response_coded = 1

# E1_3
Entrez_pheno[[6]]$Response_coded = NA
Entrez_pheno[[6]]$Response_coded[Entrez_pheno[[6]]$Response == "Responder"] = 1
Entrez_pheno[[6]]$Response_coded[Entrez_pheno[[6]]$Response == "Non_responder"] = 0

# E2
Entrez_pheno[[7]]$Response_coded = NA
Entrez_pheno[[7]]$Response_coded[Entrez_pheno[[7]]$Response == "Responder"] = 1
Entrez_pheno[[7]]$Response_coded[Entrez_pheno[[7]]$Response == "Non_responder"] = 0

# E3
Entrez_pheno[[8]]$Response_coded = NA
Entrez_pheno[[8]]$Response_coded[Entrez_pheno[[8]]$Response == "Responder"] = 1
Entrez_pheno[[8]]$Response_coded[Entrez_pheno[[8]]$Response == "Non_responder"] = 0

# E4_1
Entrez_pheno[[9]]$Response_coded = NA
Entrez_pheno[[9]]$Response_coded[Entrez_pheno[[9]]$Response == "Responder"] = 1
Entrez_pheno[[9]]$Response_coded[Entrez_pheno[[9]]$Response == "Non_responder"] = 0

# E4_2
Entrez_pheno[[10]]$Response_coded = NA
Entrez_pheno[[10]]$Response_coded[Entrez_pheno[[10]]$Response == "Responder"] = 1
Entrez_pheno[[10]]$Response_coded[Entrez_pheno[[10]]$Response == "Non_responder"] = 0

# Joining in on large phenotypic dataset
for (i in 1:length(Entrez_pheno)){
  Entrez_pheno[[i]]$Treatment_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "No_intervention", 0, 1)
  Entrez_pheno[[i]]$Endo_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "Endocrine_treatment", 1, 0)
  Entrez_pheno[[i]]$Chemo_status = ifelse(
    Entrez_pheno[[i]]$Treatment == "Chemotherapy", 1, 0)
  Entrez_pheno[[i]] = Entrez_pheno[[i]] %>% 
    dplyr::select(Sample.ID, Patient.ID, Dataset, Treatment, Treatment_status, 
                  Chemo_status, Endo_status, Response_type, Response, 
                  Response_coded, Timepoint, Timepoint_coded)
}

Pheno = rbind(Entrez_pheno[[1]], Entrez_pheno[[2]], Entrez_pheno[[3]],
              Entrez_pheno[[4]], Entrez_pheno[[5]], Entrez_pheno[[6]],
              Entrez_pheno[[7]], Entrez_pheno[[8]], Entrez_pheno[[9]],
              Entrez_pheno[[10]]); rm(i)

# Platform metadata
Pheno$Platform = as.character(Pheno$Dataset)
Pheno$Platform[Pheno$Dataset == "C1"] = "GPL17303"
Pheno$Platform[Pheno$Dataset == "C2"] = "GPL14668"
Pheno$Platform[Pheno$Dataset == "C3"] = "GPL16791"
Pheno$Platform[Pheno$Dataset == "E1_1"] = "GPL96"
Pheno$Platform[Pheno$Dataset == "E1_2"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E1_3"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E2"] = "GPL13376"
Pheno$Platform[Pheno$Dataset == "E3"] = "GPL6480"
Pheno$Platform[Pheno$Dataset == "E4_1"] = "GPL10558"
Pheno$Platform[Pheno$Dataset == "E4_2"] = "GPL10558"

Pheno$Platform_comp = as.character(Pheno$Platform)
Pheno$Platform_comp[Pheno$Dataset == "C1"] = "Ampliseq"
Pheno$Platform_comp[Pheno$Dataset == "C2"] = "Agilent"
Pheno$Platform_comp[Pheno$Dataset == "C3"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E1_1"] = "Affymetrix"
Pheno$Platform_comp[Pheno$Dataset == "E1_2"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E1_3"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E2"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E3"] = "Agilent"
Pheno$Platform_comp[Pheno$Dataset == "E4_1"] = "Illumina"
Pheno$Platform_comp[Pheno$Dataset == "E4_2"] = "Illumina"

rm(Entrez_pheno); gc()

Pheno$Response = factor(Pheno$Response, levels = c("Responder", "Non_responder"),
                        labels = c("Responder", "Non_responder"))
Pheno$Treatment = factor(Pheno$Treatment, levels = c("Chemotherapy", "Endocrine_treatment"),
                         labels = c("Chemotherapy", "Endocrine_treatment"))
Pheno$Timepoint_coded = factor(Pheno$Timepoint_coded, levels = c("T1", "T2", "T2.5", "T3"),
                               labels = c("T1", "T2", "T2.5", "T3"))
Pheno$Dataset = factor(Pheno$Dataset, levels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                                                 "E4_1", "E4_2"),
                       labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                                  "E4_1", "E4_2"))
Pheno_response = Pheno %>% dplyr::filter(is.na(Response) == FALSE)
Pheno_times = Pheno %>% dplyr::filter(Timepoint_coded == "T1" | Timepoint_coded == "T2")
Pheno_exprs = Pheno_response %>% 
  dplyr::filter(Timepoint_coded == "T1" | Timepoint_coded == "T2")

# Stratified splitting into training, validation and test #####
RNGversion("4.0.2")
set.seed(123)
tobesplit = Pheno_exprs %>%
  mutate(nr = row_number()) %>%
  dplyr::select(nr, everything()) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
train_set = tobesplit %>%
  group_by(Response, Timepoint_coded, Dataset) %>%
  sample_frac(0.7) %>%
  as.data.frame()
RNGversion("4.0.2")
set.seed(123)
validation_set = anti_join(tobesplit, train_set) %>%
  group_by(Response, Timepoint_coded, Dataset) %>%
  sample_frac(0.65) %>%
  as.data.frame()
test_set = anti_join(tobesplit, as.data.frame(rbind(train_set, validation_set)))

train_set = train_set %>% dplyr::ungroup()
validation_set = validation_set %>% dplyr::ungroup()
test_set = test_set %>% dplyr::ungroup()

train_samples = train_set$Sample.ID
validation_samples = validation_set$Sample.ID
test_samples = test_set$Sample.ID
trainval_samples = c(train_samples, validation_samples)

GEOsets_trainval = list(); GEOsets_test = list()
# For the test sets in the first loop we also add the first sample of each GEOset
# object which was dropped in lines 60-68
for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_trainval[[i]] = GEOsets[[i]][, intersect(trainval_samples, colnames(GEOsets[[i]]))]
  GEOsets_test[[i]] = GEOsets[[i]][, c(colnames(GEOsets[[i]])[1],
                                       intersect(test_samples, colnames(GEOsets[[i]])))]
}
for (i in c(1, 3, 7)){
  GEOsets_trainval[[i]] = GEOsets[[i]][, c("X", 
                                           intersect(trainval_samples, colnames(GEOsets[[i]])))]
  GEOsets_test[[i]] = GEOsets[[i]][, c("X", 
                                       intersect(test_samples, colnames(GEOsets[[i]])))]
}
names(GEOsets_test) = names(GEOsets_trainval) = names(GEOsets)
rm(i); gc()

# Probe removal (missing values) and kNN-imputation #####
# Rows that contain at least 25% of missing values are filtered out

for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_trainval[[i]] = GEOsets_trainval[[i]][rowSums(is.na(GEOsets_trainval[[i]]@assayData[["exprs"]]))/
                                length(colnames(GEOsets_trainval[[i]]@assayData[["exprs"]])) < 0.25, ]
}
rm(i)
for(i in c(1, 3, 7)){
  matrix = as.matrix(GEOsets_trainval[[i]][,2:ncol(GEOsets_trainval[[i]])])
  rownames(matrix) = GEOsets_trainval[[i]]$X
  class(matrix) = "numeric"
  matrix = matrix[rowSums(is.na(matrix))/
                    length(colnames(matrix)) < 0.25, ]
  GEOsets_trainval[[i]] = as.data.frame(matrix)
  GEOsets_trainval[[i]]$X = rownames(matrix)
  GEOsets_trainval[[i]] = GEOsets_trainval[[i]] %>% dplyr::select(X, everything())
}
rm(matrix)

# The rest of the missing values are imputed using kNN for k = 10
GEO_knn_pre1 = list(GEOsets_trainval[["C2"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_1"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_2"]]@assayData[["exprs"]],
               GEOsets_trainval[["E1_3"]]@assayData[["exprs"]],
               GEOsets_trainval[["E3"]]@assayData[["exprs"]],
               GEOsets_trainval[["E4_1"]]@assayData[["exprs"]],
               GEOsets_trainval[["E4_2"]]@assayData[["exprs"]])

GEO_knn_pre2 = list(GEOsets_trainval[["C1"]],
                    GEOsets_trainval[["C3"]],
                    GEOsets_trainval[["E2"]])

# kNN
for(i in 1:length(GEO_knn_pre1)){
  RNGversion("4.0.2")
  GEO_knn_pre1[[i]] = impute.knn(GEO_knn_pre1[[i]], k = 10, maxp = nrow(GEO_knn_pre1[[i]]),
                            rng.seed = 123)
  GEO_knn_pre1[[i]] = GEO_knn_pre1[[i]][["data"]]
}
rm(i)
for(i in 1:length(GEO_knn_pre2)){
  matrix = as.matrix(GEO_knn_pre2[[i]][,2:ncol(GEO_knn_pre2[[i]])])
  rownames(matrix) = GEO_knn_pre2[[i]]$X
  class(matrix) = "numeric"
  RNGversion("4.0.2")
  matrix = impute.knn(matrix, k = 10, maxp = nrow(matrix),
                      rng.seed = 123)[["data"]]
  GEO_knn_pre2[[i]] = as.data.frame(matrix)
  GEO_knn_pre2[[i]]$X = rownames(matrix)
  GEO_knn_pre2[[i]] = GEO_knn_pre2[[i]] %>% dplyr::select(X, everything())
}
rm(i, matrix)
names(GEO_knn_pre1) = names(GEOsets)[c(2, 4, 5, 6, 8, 9, 10)]
names(GEO_knn_pre2) = names(GEOsets)[c(1, 3, 7)]

# In order to make the test set consistent with the trainval sets we need to
# ensure the same probes are kept per study and appropriate kNN-imputation is
# performed:

for (i in c(2, 4, 5, 6, 8, 9, 10)){
  GEOsets_test[[i]] = GEOsets_test[[i]][rownames(GEOsets_trainval[[i]]), ]
}
rm(i)
for(i in c(1, 3, 7)){
  GEOsets_test[[i]] = GEOsets_test[[i]] %>%
    dplyr::filter(X %in% intersect(GEOsets_trainval[[i]]$X, GEOsets_test[[i]]$X))
}
rm(i)

# In order to impute the test samples for each study, we do the following in order:
# 1. find the 10 nearest neighbors for every gene in our imputed training and validation
#    samples of a study.
# 2. for each sample, impute the missing value in a gene using 10-NN averaging of the 
#    values in that gene's neighbours, as these were determined in step 1.
# 3. repeat for all studies

# For the external validation sets, we normalise within each study later

# The following processes are very expensive memory-wise and require a supercluster
# to run. Locally, they can be run one at a time. The loop needs excessive memory

neighbors_list_1 = list()
for (i in 1:length(GEO_knn_pre1)){
  neighbors = as.data.frame(matrix(0, ncol = 1, nrow = 10))
  colnames(neighbors) = "scaffold"
  dists = dist(GEO_knn_pre1[[i]], method = "euclidean"); gc()
  dists2 = as.matrix(dists); rm(dists); gc()
  for (j in 1:nrow(dists2)){
    gene_j_neighbors = names(sort(dists2[j,], decreasing = FALSE)[2:11]) # 2:11 is used for the top 10 neighbor genes, because the first one will always be the gene itself
    neighbors = cbind(neighbors, gene_j_neighbors)
  }
  gc()
  neighbors = neighbors[,-1]
  colnames(neighbors) = rownames(dists2)
  neighbors_list_1[[i]] = neighbors
  rm(neighbors, dists2); gc()
}
names(neighbors_list_1) = names(GEO_knn_pre1)

neighbors_list_2 = list()
for (i in 1:length(GEO_knn_pre2)){
  neighbors = as.data.frame(matrix(0, ncol = 1, nrow = 10))
  colnames(neighbors) = "scaffold"
  gem = as.matrix(GEO_knn_pre2[[i]][,-1]); rownames(gem) = GEO_knn_pre2[[i]]$X
  dists = dist(gem, method = "euclidean"); gc()
  dists2 = as.matrix(dists); rm(dists); gc()
  for (j in 1:nrow(dists2)){
    gene_j_neighbors = names(sort(dists2[j,], decreasing = FALSE)[2:11]) # 2:11 is used for the top 10 neighbor genes, because the first one will always be the gene itself
    neighbors = cbind(neighbors, gene_j_neighbors)
  }
  gc()
  neighbors = neighbors[,-1]
  colnames(neighbors) = rownames(dists2)
  neighbors_list_2[[i]] = neighbors
  rm(neighbors, dists2); gc()
}
names(neighbors_list_2) = names(GEO_knn_pre2)

# Imputing the sets that correspond to neighbors_list_1:
GEOsets_test_new = list()
for (i in c(2, 4, 5, 6, 8, 9, 10)) {
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = rownames(GEOsets_test[[i]])
  imp_map = neighbors_list_1[[names(GEOsets_test)[i]]]
  for (j in 1:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]]@assayData[["exprs"]][, j]
    for (k in 1:length(sample)) {
      if (is.na(sample[k]) == TRUE) {
        probe = names(sample[k])
        neighbors = imp_map[, probe]
        neighbor_values = sample[neighbors]
        if (sum(is.na(neighbor_values)) == 10) {
          sample[k] = mean(sample, na.rm = TRUE)
        } else {
          sample[k] = mean(neighbor_values, na.rm = TRUE)
        }
      }
    }
    scaffold = cbind(scaffold, sample)
    rm(sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  class(GEOsets_test_new[[i]]) = "numeric"
  rownames(GEOsets_test_new[[i]]) = rownames(GEOsets_test[[i]])
  rm(scaffold, probe, neighbor_values, neighbors)
}
rm(i, j, k); gc()

# Imputing the sets that correspond to neighbors_list_2:
for (i in c(1, 3, 7)) {
  scaffold = matrix(GEOsets_test[[i]]$X)
  colnames(scaffold) = "X"
  imp_map = neighbors_list_2[[names(GEOsets_test)[i]]]
  for (j in 2:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]][, j]
    names(sample) = GEOsets_test[[i]]$X
    for (k in 1:length(sample)) {
      if (is.na(sample[k]) == TRUE) {
        probe = names(sample[k])
        neighbors = imp_map[, probe]
        neighbor_values = sample[neighbors]
        if (sum(is.na(neighbor_values)) == 10) {
          sample[k] = mean(sample, na.rm = TRUE)
        } else {
          sample[k] = mean(neighbor_values, na.rm = TRUE)
        }
      }
    }
    scaffold = cbind(scaffold, sample)
    rm(sample)
  }
  GEOsets_test_new[[i]] = as.data.frame(scaffold)
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  rownames(GEOsets_test_new[[i]]) = GEOsets_test_new[[i]]$X
  GEOsets_test_new[[i]][colnames(GEOsets_test_new[[i]])[2:ncol(GEOsets_test_new[[i]])]] = 
   lapply(GEOsets_test_new[[i]][colnames(GEOsets_test_new[[i]])[2:ncol(GEOsets_test_new[[i]])]], 
         function(x) as.numeric(x))
  rm(scaffold, probe, neighbor_values, neighbors)
}
rm(i, j, k, imp_map, gene_j_neighbors); gc()
names(GEOsets_test_new) = names(GEOsets_test)

################################################################################


# The following code can be used if the test samples are to be imputed one 
# by one through an "add a sample -> impute -> remove sample" manner

# Second test set imputation approach - NOT USED IN OUR CASE #####
# Imputing the test samples one by one (for each study separately)
GEOsets_test_new = list()
for (i in c(2, 4, 5, 6, 8, 9, 10)){
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = rownames(GEOsets_test[[i]])
  impset = GEO_knn_pre1[[names(GEOsets_test)[i]]]
  for (j in 1:ncol(GEOsets_test[[i]])) {
    sample = GEOsets_test[[i]]@assayData[["exprs"]][, j]
    impmatrix = cbind(impset, sample)
    RNGversion("4.0.2")
    imp = impute.knn(impmatrix, k = 10, maxp = nrow(impmatrix),
                                   rng.seed = 123)
    imp = imp[["data"]]
    imp_sample = imp[, ncol(imp)]
    scaffold = cbind(scaffold, imp_sample)
    rm(sample, impmatrix, imp, imp_sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])
  class(GEOsets_test_new[[i]]) = "numeric"
  rm(scaffold, impset)
}
rm(i, j)

for (i in c(1, 3, 7)){
  scaffold = matrix(0, nrow = nrow(GEOsets_test[[i]]), ncol = 1)
  rownames(scaffold) = GEOsets_test[[i]]$X
  impset = as.matrix(GEO_knn_pre2[[names(GEOsets_test)[i]]][,2:ncol(GEO_knn_pre2[[names(GEOsets_test)[i]]])])
  rownames(impset) = GEO_knn_pre2[[names(GEOsets_test)[i]]]$X
  class(impset) = "numeric"
  for (j in 1:(ncol(GEOsets_test[[i]]) - 1)) {
    sample = GEOsets_test[[i]][, (j + 1)]
    impmatrix = cbind(impset, sample)
    RNGversion("4.0.2")
    imp = impute.knn(impmatrix, k = 10, maxp = nrow(impmatrix),
                     rng.seed = 123)
    imp = imp[["data"]]
    imp_sample = imp[, ncol(imp)]
    scaffold = cbind(scaffold, imp_sample)
    rm(sample, impmatrix, imp, imp_sample)
  }
  GEOsets_test_new[[i]] = scaffold[, 2:ncol(scaffold)]
  X = rownames(GEOsets_test_new[[i]])
  colnames(GEOsets_test_new[[i]]) = colnames(GEOsets_test[[i]])[2:ncol(GEOsets_test[[i]])]
  GEOsets_test_new[[i]] = as.data.frame(GEOsets_test_new[[i]])
  GEOsets_test_new[[i]]$X = X; rm(X)
  GEOsets_test_new[[i]] = GEOsets_test_new[[i]] %>% dplyr::select(X, everything())
  rm(scaffold, impset)
}

rm(i, j)
names(GEOsets_test_new) = names(GEOsets_test)

################################################################################

# Annotation: files and function ######
# ID_Map
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]

alias = org.Hs.egALIAS2EG
mapped_genes_alias = mappedkeys(alias)
alias_df = as.data.frame(alias[mapped_genes_alias])
alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
alias_df$HGNC_Official = "No"

ID_Map = rbind(official_df, alias_df) %>% distinct()
ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
  dplyr::rename(probe=Gene.Symbol) %>%
  dplyr::select(probe, EntrezGene.ID, HGNC_Official)

# Aliases
aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
Aliases = official_df %>% inner_join(aliases_for_join,
                                     by = "EntrezGene.ID") %>%
  dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
  dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
                Entrez = EntrezGene.ID) %>%
  distinct()

rm(alias, alias_df, aliases_for_join, official, mapped_genes_alias, mapped_genes_official)

# Load the annotation and genefu predictions functions
source("R_scripts/Genefu_functions.R")

# Annotation: loaded sets #####
# Annotating using annot() on the imputed data
Loaded_sets_trainval = list(GEO_knn_pre2[["C1"]], GEO_knn_pre2[["C3"]],
                            GEO_knn_pre2[["E2"]])
names(Loaded_sets_trainval) = c("C1", "C3", "E2")
Loaded_sets_test = list(GEOsets_test_new[["C1"]], GEOsets_test_new[["C3"]],
                   GEOsets_test_new[["E2"]])
names(Loaded_sets_test) = c("C1", "C3", "E2")

Loaded_exprs_Entrez_trainval = list()
Loaded_exprs_HGNC_trainval = list()
Entrez_dannot_trainval = list()
HGNC_dannot_trainval = list()

Loaded_exprs_Entrez_test = list()
Loaded_exprs_HGNC_test = list()
Entrez_dannot_test = list()
HGNC_dannot_test = list()

for (i in 1:length(Loaded_sets_trainval)){
  # Entrez-annotated data frames
  annotation = annot(expressionDF = Loaded_sets_trainval[[i]], 
                     ID_Map = ID_Map, Aliases = Aliases)
  Loaded_exprs_Entrez_trainval[[i]] = annotation$Entrez
  Loaded_exprs_Entrez_trainval[[i]] = Loaded_exprs_Entrez_trainval[[i]] %>% 
    dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
    group_by(EntrezGene.ID) %>%
    summarise_all(mean, na.rm = TRUE)
  
  Entrez_dannot_trainval[[i]] = as.data.frame(Loaded_exprs_Entrez_trainval[[i]][,"EntrezGene.ID"])
  Entrez_dannot_trainval[[i]]$probe = Entrez_dannot_trainval[[i]]$EntrezGene.ID
  Entrez_dannot_trainval[[i]] = Entrez_dannot_trainval[[i]] %>%
    dplyr::select(probe, EntrezGene.ID)
  rownames(Entrez_dannot_trainval[[i]]) = Entrez_dannot_trainval[[i]]$probe
  
  # HGNC-annotated data frames
  Loaded_exprs_HGNC_trainval[[i]] = Loaded_exprs_Entrez_trainval[[i]] %>% inner_join(official_df, by = "EntrezGene.ID") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
    dplyr::filter(is.na(Gene.Symbol) == FALSE) %>%
    dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
    group_by(Gene.Symbol) %>%
    summarise_all(mean, na.rm = TRUE)
  
  HGNC_dannot_trainval[[i]] = as.data.frame(Loaded_exprs_HGNC_trainval[[i]][,"Gene.Symbol"])
  HGNC_dannot_trainval[[i]]$probe = HGNC_dannot_trainval[[i]]$Gene.Symbol
  HGNC_dannot_trainval[[i]] = HGNC_dannot_trainval[[i]] %>%
    dplyr::select(probe, Gene.Symbol)
  rownames(HGNC_dannot_trainval[[i]]) = HGNC_dannot_trainval[[i]]$probe
}

for (i in 1:length(Loaded_sets_test)){
  # Entrez-annotated data frames
  annotation = annot(expressionDF = Loaded_sets_test[[i]], 
                     ID_Map = ID_Map, Aliases = Aliases)
  Loaded_exprs_Entrez_test[[i]] = annotation$Entrez
  Loaded_exprs_Entrez_test[[i]] = Loaded_exprs_Entrez_test[[i]] %>% 
    dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
    group_by(EntrezGene.ID) %>%
    summarise_all(mean, na.rm = TRUE)
  
  Entrez_dannot_test[[i]] = as.data.frame(Loaded_exprs_Entrez_test[[i]][,"EntrezGene.ID"])
  Entrez_dannot_test[[i]]$probe = Entrez_dannot_test[[i]]$EntrezGene.ID
  Entrez_dannot_test[[i]] = Entrez_dannot_test[[i]] %>%
    dplyr::select(probe, EntrezGene.ID)
  rownames(Entrez_dannot_test[[i]]) = Entrez_dannot_test[[i]]$probe
  
  # HGNC-annotated data frames
  Loaded_exprs_HGNC_test[[i]] = Loaded_exprs_Entrez_test[[i]] %>% inner_join(official_df, by = "EntrezGene.ID") %>%
    dplyr::select(EntrezGene.ID, Gene.Symbol, everything()) %>%
    dplyr::filter(is.na(Gene.Symbol) == FALSE) %>%
    dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
    group_by(Gene.Symbol) %>%
    summarise_all(mean, na.rm = TRUE)
  
  HGNC_dannot_test[[i]] = as.data.frame(Loaded_exprs_HGNC_test[[i]][,"Gene.Symbol"])
  HGNC_dannot_test[[i]]$probe = HGNC_dannot_test[[i]]$Gene.Symbol
  HGNC_dannot_test[[i]] = HGNC_dannot_test[[i]] %>%
    dplyr::select(probe, Gene.Symbol)
  rownames(HGNC_dannot_test[[i]]) = HGNC_dannot_test[[i]]$probe
}
rm(i, pam50, annotation); gc()
names(Loaded_exprs_Entrez_test) = names(Loaded_exprs_Entrez_trainval) = 
  names(Loaded_exprs_HGNC_test) = names(Loaded_exprs_HGNC_trainval) = 
  names(Entrez_dannot_test) = names(Entrez_dannot_trainval) = 
  names(HGNC_dannot_test) = names(HGNC_dannot_trainval) = 
  names(Loaded_sets_test)

# Annotation: sets from GEOquery (sets with whole experiments in the environment) #####
# 5 kinds of objects are created here for each set: 
# 1) exprs# is the clean expression matrix (probes matching to unique genes) - deleted later, 
# 2) exprs#_Entrez is the clean expression matrix with Entrez ID's
# 3) dannot#_Entrez is the Entrez annotation data frame for that matrix,
# 4) exprs#_HGNC is the clean expression matrix with gene symbols and
# 5) dannot#_HGNC is the HGNC annotation data frame for that second matrix

# GSE32603 - C2
exprs32603_trainval = as.data.frame(GEO_knn_pre1[["C2"]])
exprs32603_trainval$X = GEOsets_trainval[["C2"]]@featureData@data[["GENE SYMBOL"]]
exprs32603_trainval = exprs32603_trainval %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID
exprs32603_Entrez_trainval = annot(expressionDF = exprs32603_trainval, ID_Map = ID_Map,
                          Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID

dannot32603_Entrez_trainval = as.data.frame(exprs32603_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot32603_Entrez_trainval) = "EntrezGene.ID"
dannot32603_Entrez_trainval$probe = dannot32603_Entrez_trainval$EntrezGene.ID
dannot32603_Entrez_trainval = dannot32603_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot32603_Entrez_trainval) = dannot32603_Entrez_trainval$probe

# In the following chunck of code we use official_df to map Entrez ID's to official 
# Gene Symbols which are used by the Integrative Cluster algorithm. Aliases could mislead
# the algorithm in terms of ID's

exprs32603_HGNC_trainval = exprs32603_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot32603_HGNC_trainval = as.data.frame(exprs32603_HGNC_trainval[,"Gene.Symbol"])
dannot32603_HGNC_trainval$probe = dannot32603_HGNC_trainval$Gene.Symbol
dannot32603_HGNC_trainval = dannot32603_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot32603_HGNC_trainval) = dannot32603_HGNC_trainval$probe

rm(exprs32603_trainval)

# test set samples
exprs32603_test = as.data.frame(GEOsets_test_new[["C2"]])
exprs32603_test$X = GEOsets_test[["C2"]]@featureData@data[["GENE SYMBOL"]]
exprs32603_test = exprs32603_test %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID
exprs32603_Entrez_test = annot(expressionDF = exprs32603_test, ID_Map = ID_Map,
                               Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE) # get one row for each gene which has the mean of all probes matching to that ID

dannot32603_Entrez_test = as.data.frame(exprs32603_Entrez_test[,"EntrezGene.ID"])
colnames(dannot32603_Entrez_test) = "EntrezGene.ID"
dannot32603_Entrez_test$probe = dannot32603_Entrez_test$EntrezGene.ID
dannot32603_Entrez_test = dannot32603_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot32603_Entrez_test) = dannot32603_Entrez_test$probe

# In the following chunck of code we use official_df to map Entrez ID's to official 
# Gene Symbols which are used by the Integrative Cluster algorithm. Aliases could mislead
# the algorithm in terms of ID's

exprs32603_HGNC_test = exprs32603_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot32603_HGNC_test = as.data.frame(exprs32603_HGNC_test[,"Gene.Symbol"])
dannot32603_HGNC_test$probe = dannot32603_HGNC_test$Gene.Symbol
dannot32603_HGNC_test = dannot32603_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot32603_HGNC_test) = dannot32603_HGNC_test$probe

rm(exprs32603_test)

# GSE20181 - E1_1
feature_trainval = GEOsets_trainval[["E1_1"]]@featureData@data %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
exprs20181_trainval = as.data.frame(GEO_knn_pre1[["E1_1"]])[feature_trainval$ID,]
exprs20181_trainval$ID = rownames(exprs20181_trainval)
exprs20181_Entrez_trainval = exprs20181_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  dplyr::filter(is.na(ENTREZ_GENE_ID) == FALSE) %>% # just in case
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = ENTREZ_GENE_ID)

dannot20181_Entrez_trainval = as.data.frame(exprs20181_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot20181_Entrez_trainval) = "EntrezGene.ID"
dannot20181_Entrez_trainval$probe = dannot20181_Entrez_trainval$EntrezGene.ID
dannot20181_Entrez_trainval = dannot20181_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot20181_Entrez_trainval) = dannot20181_Entrez_trainval$probe

exprs20181_HGNC_trainval = exprs20181_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot20181_HGNC_trainval = as.data.frame(exprs20181_HGNC_trainval[,"Gene.Symbol"])
dannot20181_HGNC_trainval$probe = dannot20181_HGNC_trainval$Gene.Symbol
dannot20181_HGNC_trainval = dannot20181_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot20181_HGNC_trainval) = dannot20181_HGNC_trainval$probe

rm(exprs20181_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_1"]]@featureData@data %>%
  dplyr::select(ID, ENTREZ_GENE_ID) %>%
  dplyr::filter(!grepl("///", ENTREZ_GENE_ID)) %>%
  dplyr::filter(nchar(ENTREZ_GENE_ID)>0)
exprs20181_test = as.data.frame(GEOsets_test_new[["E1_1"]])[feature_test$ID,]
exprs20181_test$ID = rownames(exprs20181_test)
exprs20181_Entrez_test = exprs20181_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(ENTREZ_GENE_ID, everything()) %>%
  dplyr::filter(is.na(ENTREZ_GENE_ID) == FALSE) %>% # just in case
  group_by(ENTREZ_GENE_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = ENTREZ_GENE_ID)

dannot20181_Entrez_test = as.data.frame(exprs20181_Entrez_test[,"EntrezGene.ID"])
colnames(dannot20181_Entrez_test) = "EntrezGene.ID"
dannot20181_Entrez_test$probe = dannot20181_Entrez_test$EntrezGene.ID
dannot20181_Entrez_test = dannot20181_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot20181_Entrez_test) = dannot20181_Entrez_test$probe

exprs20181_HGNC_test = exprs20181_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot20181_HGNC_test = as.data.frame(exprs20181_HGNC_test[,"Gene.Symbol"])
dannot20181_HGNC_test$probe = dannot20181_HGNC_test$Gene.Symbol
dannot20181_HGNC_test = dannot20181_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot20181_HGNC_test) = dannot20181_HGNC_test$probe

rm(exprs20181_test, feature_test)

# GSE55374 - E1_2
feature_trainval = GEOsets_trainval[["E1_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs55374_trainval = as.data.frame(GEO_knn_pre1[["E1_2"]])[feature_trainval$ID,]
exprs55374_trainval$ID = rownames(exprs55374_trainval)
exprs55374_Entrez_trainval = exprs55374_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs55374_Entrez_trainval$EntrezGene.ID = as.character(exprs55374_Entrez_trainval$EntrezGene.ID)

dannot55374_Entrez_trainval = as.data.frame(exprs55374_Entrez_trainval[,"EntrezGene.ID"])
dannot55374_Entrez_trainval$probe = dannot55374_Entrez_trainval$EntrezGene.ID
dannot55374_Entrez_trainval = dannot55374_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot55374_Entrez_trainval) = dannot55374_Entrez_trainval$probe

exprs55374_HGNC_trainval = exprs55374_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot55374_HGNC_trainval = as.data.frame(exprs55374_HGNC_trainval[,"Gene.Symbol"])
dannot55374_HGNC_trainval$probe = dannot55374_HGNC_trainval$Gene.Symbol
dannot55374_HGNC_trainval = dannot55374_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot55374_HGNC_trainval) = dannot55374_HGNC_trainval$probe

rm(exprs55374_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs55374_test = as.data.frame(GEOsets_test_new[["E1_2"]])[feature_test$ID,]
exprs55374_test$ID = rownames(exprs55374_test)
exprs55374_Entrez_test = exprs55374_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs55374_Entrez_test$EntrezGene.ID = as.character(exprs55374_Entrez_test$EntrezGene.ID)

dannot55374_Entrez_test = as.data.frame(exprs55374_Entrez_test[,"EntrezGene.ID"])
dannot55374_Entrez_test$probe = dannot55374_Entrez_test$EntrezGene.ID
dannot55374_Entrez_test = dannot55374_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot55374_Entrez_test) = dannot55374_Entrez_test$probe

exprs55374_HGNC_test = exprs55374_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot55374_HGNC_test = as.data.frame(exprs55374_HGNC_test[,"Gene.Symbol"])
dannot55374_HGNC_test$probe = dannot55374_HGNC_test$Gene.Symbol
dannot55374_HGNC_test = dannot55374_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot55374_HGNC_test) = dannot55374_HGNC_test$probe

rm(exprs55374_test, feature_test)

# GSE59515 - E1_3
feature_trainval = GEOsets_trainval[["E1_3"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs59515_trainval = as.data.frame(GEO_knn_pre1[["E1_3"]])[feature_trainval$ID,]
exprs59515_trainval$ID = rownames(exprs59515_trainval)
exprs59515_Entrez_trainval = exprs59515_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs59515_Entrez_trainval$EntrezGene.ID = as.character(exprs59515_Entrez_trainval$EntrezGene.ID)

dannot59515_Entrez_trainval = as.data.frame(exprs59515_Entrez_trainval[,"EntrezGene.ID"])
dannot59515_Entrez_trainval$probe = dannot59515_Entrez_trainval$EntrezGene.ID
dannot59515_Entrez_trainval = dannot59515_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot59515_Entrez_trainval) = dannot59515_Entrez_trainval$probe

exprs59515_HGNC_trainval = exprs59515_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot59515_HGNC_trainval = as.data.frame(exprs59515_HGNC_trainval[,"Gene.Symbol"])
dannot59515_HGNC_trainval$probe = dannot59515_HGNC_trainval$Gene.Symbol
dannot59515_HGNC_trainval = dannot59515_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot59515_HGNC_trainval) = dannot59515_HGNC_trainval$probe

rm(exprs59515_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E1_3"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs59515_test = as.data.frame(GEOsets_test_new[["E1_3"]])[feature_test$ID,]
exprs59515_test$ID = rownames(exprs59515_test)
exprs59515_Entrez_test = exprs59515_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs59515_Entrez_test$EntrezGene.ID = as.character(exprs59515_Entrez_test$EntrezGene.ID)

dannot59515_Entrez_test = as.data.frame(exprs59515_Entrez_test[,"EntrezGene.ID"])
dannot59515_Entrez_test$probe = dannot59515_Entrez_test$EntrezGene.ID
dannot59515_Entrez_test = dannot59515_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot59515_Entrez_test) = dannot59515_Entrez_test$probe

exprs59515_HGNC_test = exprs59515_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot59515_HGNC_test = as.data.frame(exprs59515_HGNC_test[,"Gene.Symbol"])
dannot59515_HGNC_test$probe = dannot59515_HGNC_test$Gene.Symbol
dannot59515_HGNC_test = dannot59515_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot59515_HGNC_test) = dannot59515_HGNC_test$probe

rm(exprs59515_test, feature_test)

# GSE87411 - E3
exprs87411_trainval = as.data.frame(GEO_knn_pre1[["E3"]])
exprs87411_trainval$X = GEOsets_trainval[["E3"]]@featureData@data[["GENE_SYMBOL"]]
exprs87411_trainval = exprs87411_trainval %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE)
exprs87411_Entrez_trainval = annot(expressionDF = exprs87411_trainval, ID_Map = ID_Map,
                          Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_Entrez_trainval = as.data.frame(exprs87411_Entrez_trainval[,"EntrezGene.ID"])
colnames(dannot87411_Entrez_trainval) = "EntrezGene.ID"
dannot87411_Entrez_trainval$probe = dannot87411_Entrez_trainval$EntrezGene.ID
dannot87411_Entrez_trainval = dannot87411_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot87411_Entrez_trainval) = dannot87411_Entrez_trainval$probe

exprs87411_HGNC_trainval = exprs87411_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_HGNC_trainval = as.data.frame(exprs87411_HGNC_trainval[,"Gene.Symbol"])
dannot87411_HGNC_trainval$probe = dannot87411_HGNC_trainval$Gene.Symbol
dannot87411_HGNC_trainval = dannot87411_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot87411_HGNC_trainval) = dannot87411_HGNC_trainval$probe

rm(exprs87411_trainval)

# test set samples
exprs87411_test = as.data.frame(GEOsets_test_new[["E3"]])
exprs87411_test$X = GEOsets_test[["E3"]]@featureData@data[["GENE_SYMBOL"]]
exprs87411_test = exprs87411_test %>% dplyr::select(X, everything()) %>%
  dplyr::filter(nchar(X)>0) %>%
  dplyr::filter(!grepl("///", X)) %>% # in case there are multiple genes in a probe
  dplyr::filter(!grepl(",", X)) %>% # in case there are multiple genes in a probe
  group_by(X) %>%
  summarise_all(mean, na.rm = TRUE)
exprs87411_Entrez_test = annot(expressionDF = exprs87411_test, ID_Map = ID_Map,
                               Aliases = Aliases)$Entrez %>%
  dplyr::filter(is.na(EntrezGene.ID) == FALSE) %>%
  group_by(EntrezGene.ID) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_Entrez_test = as.data.frame(exprs87411_Entrez_test[,"EntrezGene.ID"])
colnames(dannot87411_Entrez_test) = "EntrezGene.ID"
dannot87411_Entrez_test$probe = dannot87411_Entrez_test$EntrezGene.ID
dannot87411_Entrez_test = dannot87411_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot87411_Entrez_test) = dannot87411_Entrez_test$probe

exprs87411_HGNC_test = exprs87411_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot87411_HGNC_test = as.data.frame(exprs87411_HGNC_test[,"Gene.Symbol"])
dannot87411_HGNC_test$probe = dannot87411_HGNC_test$Gene.Symbol
dannot87411_HGNC_test = dannot87411_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot87411_HGNC_test) = dannot87411_HGNC_test$probe

rm(exprs87411_test)

# GSE105777 - E4_1
feature_trainval = GEOsets_trainval[["E4_1"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs105777_trainval = as.data.frame(GEO_knn_pre1[["E4_1"]])[feature_trainval$ID,]
exprs105777_trainval$ID = rownames(exprs105777_trainval)
exprs105777_Entrez_trainval = exprs105777_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs105777_Entrez_trainval$EntrezGene.ID = as.character(exprs105777_Entrez_trainval$EntrezGene.ID)

dannot105777_Entrez_trainval = as.data.frame(exprs105777_Entrez_trainval[,"EntrezGene.ID"])
dannot105777_Entrez_trainval$probe = dannot105777_Entrez_trainval$EntrezGene.ID
dannot105777_Entrez_trainval = dannot105777_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot105777_Entrez_trainval) = dannot105777_Entrez_trainval$probe

exprs105777_HGNC_trainval = exprs105777_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot105777_HGNC_trainval = as.data.frame(exprs105777_HGNC_trainval[,"Gene.Symbol"])
dannot105777_HGNC_trainval$probe = dannot105777_HGNC_trainval$Gene.Symbol
dannot105777_HGNC_trainval = dannot105777_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot105777_HGNC_trainval) = dannot105777_HGNC_trainval$probe

rm(exprs105777_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E4_1"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs105777_test = as.data.frame(GEOsets_test_new[["E4_1"]])[feature_test$ID,]
exprs105777_test$ID = rownames(exprs105777_test)
exprs105777_Entrez_test = exprs105777_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs105777_Entrez_test$EntrezGene.ID = as.character(exprs105777_Entrez_test$EntrezGene.ID)

dannot105777_Entrez_test = as.data.frame(exprs105777_Entrez_test[,"EntrezGene.ID"])
dannot105777_Entrez_test$probe = dannot105777_Entrez_test$EntrezGene.ID
dannot105777_Entrez_test = dannot105777_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot105777_Entrez_test) = dannot105777_Entrez_test$probe

exprs105777_HGNC_test = exprs105777_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot105777_HGNC_test = as.data.frame(exprs105777_HGNC_test[,"Gene.Symbol"])
dannot105777_HGNC_test$probe = dannot105777_HGNC_test$Gene.Symbol
dannot105777_HGNC_test = dannot105777_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot105777_HGNC_test) = dannot105777_HGNC_test$probe

rm(exprs105777_test, feature_test)

# GSE126870 - E4_2
feature_trainval = GEOsets_trainval[["E4_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs126870_trainval = as.data.frame(GEO_knn_pre1[["E4_2"]])[feature_trainval$ID,]
exprs126870_trainval$ID = rownames(exprs126870_trainval)
exprs126870_Entrez_trainval = exprs126870_trainval %>%
  inner_join(feature_trainval) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs126870_Entrez_trainval$EntrezGene.ID = as.character(exprs126870_Entrez_trainval$EntrezGene.ID)

dannot126870_Entrez_trainval = as.data.frame(exprs126870_Entrez_trainval[,"EntrezGene.ID"])
dannot126870_Entrez_trainval$probe = dannot126870_Entrez_trainval$EntrezGene.ID
dannot126870_Entrez_trainval = dannot126870_Entrez_trainval %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot126870_Entrez_trainval) = dannot126870_Entrez_trainval$probe

exprs126870_HGNC_trainval = exprs126870_Entrez_trainval %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot126870_HGNC_trainval = as.data.frame(exprs126870_HGNC_trainval[,"Gene.Symbol"])
dannot126870_HGNC_trainval$probe = dannot126870_HGNC_trainval$Gene.Symbol
dannot126870_HGNC_trainval = dannot126870_HGNC_trainval %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot126870_HGNC_trainval) = dannot126870_HGNC_trainval$probe

rm(exprs126870_trainval, feature_trainval)

# test set samples
feature_test = GEOsets_test[["E4_2"]]@featureData@data %>%
  dplyr::select(ID, Entrez_Gene_ID) %>%
  dplyr::filter(!grepl("///", Entrez_Gene_ID)) %>%
  dplyr::filter(!grepl(",", Entrez_Gene_ID)) %>% # in case there are multiple genes in a probe
  dplyr::filter(nchar(Entrez_Gene_ID)>0)
exprs126870_test = as.data.frame(GEOsets_test_new[["E4_2"]])[feature_test$ID,]
exprs126870_test$ID = rownames(exprs126870_test)
exprs126870_Entrez_test = exprs126870_test %>%
  inner_join(feature_test) %>%
  dplyr::select(-ID) %>%
  dplyr::select(Entrez_Gene_ID, everything()) %>%
  dplyr::filter(is.na(Entrez_Gene_ID) == FALSE) %>% # just in case
  group_by(Entrez_Gene_ID) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::rename(EntrezGene.ID = Entrez_Gene_ID)
exprs126870_Entrez_test$EntrezGene.ID = as.character(exprs126870_Entrez_test$EntrezGene.ID)

dannot126870_Entrez_test = as.data.frame(exprs126870_Entrez_test[,"EntrezGene.ID"])
dannot126870_Entrez_test$probe = dannot126870_Entrez_test$EntrezGene.ID
dannot126870_Entrez_test = dannot126870_Entrez_test %>%
  dplyr::select(probe, EntrezGene.ID)
rownames(dannot126870_Entrez_test) = dannot126870_Entrez_test$probe

exprs126870_HGNC_test = exprs126870_Entrez_test %>% inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(-EntrezGene.ID, -HGNC_Official) %>%
  dplyr::select(Gene.Symbol, everything()) %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean, na.rm = TRUE)

dannot126870_HGNC_test = as.data.frame(exprs126870_HGNC_test[,"Gene.Symbol"])
dannot126870_HGNC_test$probe = dannot126870_HGNC_test$Gene.Symbol
dannot126870_HGNC_test = dannot126870_HGNC_test %>%
  dplyr::select(probe, Gene.Symbol)
rownames(dannot126870_HGNC_test) = dannot126870_HGNC_test$probe

rm(exprs126870_test, feature_test)

# Joining them in a list of Entrez-annotated sets and a list of HGNC-annotated sets
Entrez_exprs_trainval = list(Loaded_exprs_Entrez_trainval[["C1"]], exprs32603_Entrez_trainval, 
                             Loaded_exprs_Entrez_trainval[["C3"]], exprs20181_Entrez_trainval, 
                             exprs55374_Entrez_trainval, exprs59515_Entrez_trainval, 
                    Loaded_exprs_Entrez_trainval[["E2"]], exprs87411_Entrez_trainval, 
                    exprs105777_Entrez_trainval, exprs126870_Entrez_trainval)
names(Entrez_exprs_trainval) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                        "E4_1", "E4_2")

Entrez_exprs_test = list(Loaded_exprs_Entrez_test[["C1"]], exprs32603_Entrez_test, 
                         Loaded_exprs_Entrez_test[["C3"]], exprs20181_Entrez_test, 
                         exprs55374_Entrez_test, exprs59515_Entrez_test, 
                         Loaded_exprs_Entrez_test[["E2"]], exprs87411_Entrez_test, 
                         exprs105777_Entrez_test, exprs126870_Entrez_test)
names(Entrez_exprs_test) = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3", "E2", "E3",
                             "E4_1", "E4_2")

HGNC_exprs_trainval = list(Loaded_exprs_HGNC_trainval[["C1"]], exprs32603_HGNC_trainval, 
                           Loaded_exprs_HGNC_trainval[["C3"]], exprs20181_HGNC_trainval, 
                           exprs55374_HGNC_trainval, exprs59515_HGNC_trainval, 
                  Loaded_exprs_HGNC_trainval[["E2"]], exprs87411_HGNC_trainval, 
                  exprs105777_HGNC_trainval, exprs126870_HGNC_trainval)
names(HGNC_exprs_trainval) = names(Entrez_exprs_trainval)

HGNC_exprs_test = list(Loaded_exprs_HGNC_test[["C1"]], exprs32603_HGNC_test, 
                       Loaded_exprs_HGNC_test[["C3"]], exprs20181_HGNC_test, 
                       exprs55374_HGNC_test, exprs59515_HGNC_test, 
                       Loaded_exprs_HGNC_test[["E2"]], exprs87411_HGNC_test, 
                       exprs105777_HGNC_test, exprs126870_HGNC_test)
names(HGNC_exprs_test) = names(Entrez_exprs_test)

dannot_Entrez_trainval = list(Entrez_dannot_trainval[["C1"]], dannot32603_Entrez_trainval, 
                              Entrez_dannot_trainval[["C3"]], dannot20181_Entrez_trainval, 
                              dannot55374_Entrez_trainval, dannot59515_Entrez_trainval, 
                     Entrez_dannot_trainval[["E2"]], dannot87411_Entrez_trainval, 
                     dannot105777_Entrez_trainval, dannot126870_Entrez_trainval)
names(dannot_Entrez_trainval) = names(Entrez_exprs_trainval)

dannot_Entrez_test = list(Entrez_dannot_test[["C1"]], dannot32603_Entrez_test, 
                          Entrez_dannot_test[["C3"]], dannot20181_Entrez_test, 
                          dannot55374_Entrez_test, dannot59515_Entrez_test, 
                          Entrez_dannot_test[["E2"]], dannot87411_Entrez_test, 
                          dannot105777_Entrez_test, dannot126870_Entrez_test)
names(dannot_Entrez_test) = names(Entrez_exprs_test)

dannot_HGNC_trainval = list(HGNC_dannot_trainval[["C1"]], dannot32603_HGNC_trainval, 
                            HGNC_dannot_trainval[["C3"]], dannot20181_HGNC_trainval, 
                            dannot55374_HGNC_trainval, dannot59515_HGNC_trainval, 
                   HGNC_dannot_trainval[["E2"]], dannot87411_HGNC_trainval, 
                   dannot105777_HGNC_trainval, dannot126870_HGNC_trainval)
names(dannot_HGNC_trainval) = names(Entrez_exprs_trainval)

dannot_HGNC_test = list(HGNC_dannot_test[["C1"]], dannot32603_HGNC_test, 
                        HGNC_dannot_test[["C3"]], dannot20181_HGNC_test, 
                        dannot55374_HGNC_test, dannot59515_HGNC_test, 
                        HGNC_dannot_test[["E2"]], dannot87411_HGNC_test, 
                        dannot105777_HGNC_test, dannot126870_HGNC_test)
names(dannot_HGNC_test) = names(Entrez_exprs_test)

rm(exprs20181_Entrez_trainval, exprs32603_Entrez_trainval, exprs55374_Entrez_trainval,
   exprs59515_Entrez_trainval, exprs87411_Entrez_trainval, Loaded_exprs_Entrez_trainval, 
   exprs20181_HGNC_trainval, exprs32603_HGNC_trainval, exprs55374_HGNC_trainval, 
   exprs59515_HGNC_trainval, exprs87411_HGNC_trainval, Loaded_exprs_HGNC_trainval, 
   Loaded_sets_trainval, Entrez_dannot_trainval, HGNC_dannot_trainval, 
   dannot20181_HGNC_trainval, dannot32603_HGNC_trainval, dannot55374_HGNC_trainval, 
   dannot59515_HGNC_trainval, dannot87411_HGNC_trainval, dannot20181_Entrez_trainval, 
   dannot32603_Entrez_trainval, dannot55374_Entrez_trainval, dannot59515_Entrez_trainval, 
   dannot87411_Entrez_trainval, pam50, exprs105777_Entrez_trainval, exprs105777_HGNC_trainval, 
   exprs126870_Entrez_trainval, exprs126870_HGNC_trainval, dannot105777_Entrez_trainval, 
   dannot126870_Entrez_trainval, dannot105777_HGNC_trainval, dannot126870_HGNC_trainval,
   exprs20181_Entrez_test, exprs32603_Entrez_test, exprs55374_Entrez_test,
   exprs59515_Entrez_test, exprs87411_Entrez_test, Loaded_exprs_Entrez_test, 
   exprs20181_HGNC_test, exprs32603_HGNC_test, exprs55374_HGNC_test, 
   exprs59515_HGNC_test, exprs87411_HGNC_test, Loaded_exprs_HGNC_test, 
   Loaded_sets_test, Entrez_dannot_test, HGNC_dannot_test, 
   dannot20181_HGNC_test, dannot32603_HGNC_test, dannot55374_HGNC_test, 
   dannot59515_HGNC_test, dannot87411_HGNC_test, dannot20181_Entrez_test, 
   dannot32603_Entrez_test, dannot55374_Entrez_test, dannot59515_Entrez_test, 
   dannot87411_Entrez_test, exprs105777_Entrez_test, exprs105777_HGNC_test, 
   exprs126870_Entrez_test, exprs126870_HGNC_test, dannot105777_Entrez_test, 
   dannot126870_Entrez_test, dannot105777_HGNC_test, dannot126870_HGNC_test)
gc()

# Genefu phenotypic annotation #####

Genefu_sets_trainval = list()
for (i in 1:length(Entrez_exprs_trainval)){
  genefu_set = Genefu_predictions_2(expressionDF = Entrez_exprs_trainval[[i]],
                                  HGNC = HGNC_exprs_trainval[[i]],
                                  dannot_Entrez = dannot_Entrez_trainval[[i]],
                                  dannot_HGNC = dannot_HGNC_trainval[[i]])
  Genefu_sets_trainval[[i]] = genefu_set
}

rm(claudinLowData, genefu_set, pam50.robust, scmod1.robust, sig.gene70, i); gc()
names(Genefu_sets_trainval) = names(GEOsets)

Genefu_sets_test = list()
for (i in 1:length(Entrez_exprs_test)){
  genefu_set = Genefu_predictions_2(expressionDF = Entrez_exprs_test[[i]],
                                  HGNC = HGNC_exprs_test[[i]],
                                  dannot_Entrez = dannot_Entrez_test[[i]],
                                  dannot_HGNC = dannot_HGNC_test[[i]])
  Genefu_sets_test[[i]] = genefu_set
}

rm(claudinLowData, genefu_set, pam50.robust, scmod1.robust, sig.gene70, i); gc()
names(Genefu_sets_test) = names(GEOsets)

Genefu_megaset_trainval = rbind(Genefu_sets_trainval[["C1"]], Genefu_sets_trainval[["C2"]], 
                                Genefu_sets_trainval[["C3"]], Genefu_sets_trainval[["E1_1"]], 
                                Genefu_sets_trainval[["E1_2"]], Genefu_sets_trainval[["E1_3"]],
                                Genefu_sets_trainval[["E2"]], Genefu_sets_trainval[["E3"]],
                                Genefu_sets_trainval[["E4_1"]], Genefu_sets_trainval[["E4_2"]])

Genefu_megaset_test = rbind(Genefu_sets_test[["C1"]], Genefu_sets_test[["C2"]], 
                            Genefu_sets_test[["C3"]], Genefu_sets_test[["E1_1"]], 
                            Genefu_sets_test[["E1_2"]], Genefu_sets_test[["E1_3"]],
                            Genefu_sets_test[["E2"]], Genefu_sets_test[["E3"]],
                            Genefu_sets_test[["E4_1"]], Genefu_sets_test[["E4_2"]])

gc()

Pheno_trainval = rbind(train_set, validation_set) %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")
Pheno_train = train_set %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")
Pheno_val = validation_set %>% 
  inner_join(Genefu_megaset_trainval, by = "Sample.ID")

# Manually append test set with the dropped samples prior to joining:
rogues = setdiff(Genefu_megaset_test$Sample.ID, test_set$Sample.ID)
one = c("5105", rogues[1], "C2.1001", "C2", "Chemotherapy", "1", "1", "0",
        "RCB, cutoff: >1", "Non_responder", "0", "Pre-treatment", "T1", "GPL14668",
        "Agilent")
two = c("5145", rogues[2], "E1.1.10", "E1_1", "Endocrine_treatment", "1", "0", "1",
        "US_tumour_volume", "Responder", "1", "Pre-treatment", "T1", "GPL96",
        "Affymetrix")
three = c("6145", rogues[3], "E1.2.182", "E1_2", "Endocrine_treatment", "1", "0", "1",
          "US_tumour_volume", "Responder", "1", "Pre-treatment", "T1", "GPL10558",
          "Illumina")
four = c("4145", rogues[4], "E1.3.150", "E1_3", "Endocrine_treatment", "1", "0", "1",
         "US_tumour_volume", "Non_responder", "0", "Pre-treatment", "T1", "GPL10558",
         "Illumina")
five = c("4195", rogues[5], "E3.16723", "E3", "Endocrine_treatment", "1", "0", "1",
         "pCR (+Ki67)", "Non_responder", "0", "Pre-treatment", "T1", "GPL6480",
         "Agilent")
six = c("4159", rogues[6], "E4.1", "E4_1", "Endocrine_treatment", "1", "0", "1",
        "Ki67 60% expr. decr.", "Responder", "1", "Pre-treatment", "T1", "GPL10558",
        "Illumina")
seven = c("4165", rogues[7], "E4.174", "E4_2", "Endocrine_treatment", "1", "0", "1",
          "Ki67 60% expr. decr.", "Non_responder", "0", "Pre-treatment", "T1", "GPL10558",
          "Illumina")

test_set = rbind(test_set, one, two, three, four, five, six, seven)
Pheno = rbind(Pheno, one, two, three, four, five, six, seven)
Pheno_test = test_set %>%
  inner_join(Genefu_megaset_test, by = "Sample.ID")
Pheno_full = rbind(Pheno_trainval, Pheno_test)
Pheno_full = Pheno_full[order(Pheno_full$Dataset),]
rm(Genefu_megaset_trainval, Genefu_megaset_test, Genefu_sets_trainval, Genefu_sets_test,
   dannot_Entrez_trainval, dannot_Entrez_test, dannot_HGNC_trainval,
   dannot_HGNC_test, HGNC_exprs_trainval, HGNC_exprs_test, one, two, three, four, five,
   six, seven, rogues); gc()

Pairwise = inner_join(Pheno[which(Pheno$Timepoint_coded == "T1"),],
                      Pheno[which(Pheno$Timepoint_coded == "T2"),], 
                      by = c("Patient.ID", "Treatment", "Treatment_status",
                             "Chemo_status", "Endo_status", "Response_type", "Response",
                             "Response_coded", "Platform", "Platform_comp")) %>%
  dplyr::select(-Timepoint.y, -Timepoint_coded.y, -Timepoint.x, -Timepoint_coded.x)

# Overview of samples, patients and pre-/on-pairs
rbind(table(Pheno$Dataset), cbind(length(unique(Pheno$Patient.ID[Pheno$Dataset == "C1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "C2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "C3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E1_3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E2"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E3"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E4_1"])), 
                                  length(unique(Pheno$Patient.ID[Pheno$Dataset == "E4_2"]))), 
      table(Pairwise$Dataset.x), table(Pairwise$Dataset.y))

# Taking into account the fact that the last two datasets have some overlap, we have:
# 1610 (out of original 1858) samples, 824 patients and 545 pre-/on-pairs

phenowb = createWorkbook()
addWorksheet(phenowb, "Full")
writeData(phenowb, "Full", Pheno_full)
addWorksheet(phenowb, "Train_val")
writeData(phenowb, "Train_val", Pheno_trainval)
addWorksheet(phenowb, "Train")
writeData(phenowb, "Train", Pheno_train)
addWorksheet(phenowb, "Validation")
writeData(phenowb, "Validation", Pheno_val)
addWorksheet(phenowb, "Test")
writeData(phenowb, "Test", Pheno_test)
saveWorkbook(phenowb, "data/Output sets/Pheno.xlsx", overwrite = TRUE)


detach("package:GEOquery", unload = TRUE)
detach("package:readr", unload = TRUE)
library(readr)
# QC plots #####
# Preparation

library(reshape2)
library(pheatmap)
library(ggplot2)

exprs = inner_join(Entrez_exprs_trainval[[1]], Entrez_exprs_trainval[[2]], 
                   by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[3]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[4]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[5]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[6]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[7]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[8]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[9]], by = "EntrezGene.ID") %>%
  inner_join(Entrez_exprs_trainval[[10]], by = "EntrezGene.ID")

exprs = na.omit(exprs)
Pheno_exprs = Pheno_trainval
check_exprs = as.matrix(exprs[, 2:ncol(exprs)])
check_exprs = check_exprs[, Pheno_exprs$Sample.ID]
rownames(check_exprs) = exprs$EntrezGene.ID

plot_matrix = as.matrix(exprs[,2:ncol(exprs)])
plot_matrix = plot_matrix[, Pheno_exprs$Sample.ID]
rownames(plot_matrix) = exprs$EntrezGene.ID
mds = plotMDS(plot_matrix)
pca = data.frame(cbind(mds$x, mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                       Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                       as.character(Pheno_exprs$Timepoint_coded),
                       Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                       as.character(Pheno_exprs$IC10), 
                       as.character(Pheno_exprs$Mammaprint_risk), 
                       as.character(Pheno_exprs$rorS_risk)))
pca$Mammaprint = pca$X11
pca$Mammaprint[which(pca$X11=="1")] = "Risk"
pca$Mammaprint[which(pca$X11!="1")] = "No risk"
pca = pca %>% dplyr::select(-X11)
colnames(pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                  "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  pca[,i] = factor(pca[,i], levels = unique(pca[,i]), labels = unique(pca[,i]))
}

class(pca$X1) = "numeric"
class(pca$X2) = "numeric"

# Multidimensional scaling plots #####
gpca = ggplot(pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS22 (", round(100*mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
gpca
ggsave(filename = "MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# RBE
RBE_design = model.matrix(~0 + Pheno_exprs$Response + 
                            Pheno_exprs$Timepoint_coded +
                            Pheno_exprs$pam50)[,c(1:3, 6:9)]
colnames(RBE_design) = c("Non_responder", "Responder", "T2",
                         "HER2", "LumB", "LumA", "Normal")
rownames(RBE_design) = Pheno_exprs$Sample.ID
cm_RBE = makeContrasts(RespvsNonresp = Responder - Non_responder,
                       levels = RBE_design)
seq_batch = Pheno_exprs$Dataset
RBE_exprs = removeBatchEffect(check_exprs, batch = seq_batch, design = RBE_design)

RBE_mds = plotMDS(RBE_exprs)
RBE_pca = data.frame(cbind(RBE_mds$x, RBE_mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                           Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                           as.character(Pheno_exprs$Timepoint_coded),
                           Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                           as.character(Pheno_exprs$IC10), 
                           as.character(Pheno_exprs$Mammaprint_risk), 
                           as.character(Pheno_exprs$rorS_risk)))
RBE_pca$Mammaprint = RBE_pca$X11
RBE_pca$Mammaprint[which(RBE_pca$X11=="1")] = "Risk"
RBE_pca$Mammaprint[which(RBE_pca$X11!="1")] = "No risk"
RBE_pca = RBE_pca %>% dplyr::select(-X11)
colnames(RBE_pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                      "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  RBE_pca[,i] = factor(RBE_pca[,i], levels = unique(RBE_pca[,i]), labels = unique(RBE_pca[,i]))
}

class(RBE_pca$X1) = "numeric"
class(RBE_pca$X2) = "numeric"

rbe_gpca = ggplot(RBE_pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "RBE - Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*RBE_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS2 (", round(100*RBE_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
rbe_gpca
ggsave(filename = "RBE_MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# z-score normalisation ( https:://www.biostars.org/p/283083/ )
z = list()
means = list()
sd = list()

for(i in 1:length(Entrez_exprs_trainval)){
  df = as.data.frame(Entrez_exprs_trainval[[i]][,2:ncol(Entrez_exprs_trainval[[i]])])
  t = as.data.frame(t(df))
  means[[i]] = sapply(t, function(t) mean(t, na.rm = T))
  sd[[i]] = sapply(t, function(t) sd(t, na.rm = T))
  z_t = sapply(t, function(t) (t-mean(t, na.rm = T))/sd(t, na.rm = T))
  z[[i]] = as.matrix(t(z_t))
  rownames(z[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
  colnames(z[[i]]) = colnames(df)
  z[[i]] = as.data.frame(z[[i]])
  z[[i]]$EntrezGene.ID = Entrez_exprs_trainval[[i]]$EntrezGene.ID
  rm(t, z_t, df)
}

for (i in 1:length(means)) names(means[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
for (i in 1:length(sd)) names(sd[[i]]) = Entrez_exprs_trainval[[i]]$EntrezGene.ID
names(z) = names(means) = names(sd) = names(Entrez_exprs_trainval)

# To apply normalisation consistently on the test set we will use the saved 
# means and standard deviations in the means and sd lists that we created:

z_test = list()
for (i in 1:length(Entrez_exprs_test)){
  df = as.data.frame(Entrez_exprs_test[[i]][,2:ncol(Entrez_exprs_test[[i]])])
  t = as.data.frame(t(df))
  z_t = t
  for (j in 1:nrow(t)){
    for (k in 1:ncol(t)){
      z_t[j, k] = (z_t[j, k] - means[[i]][k])/sd[[i]][k]
    }
  }
  z_test[[i]] = as.matrix(t(z_t))
  rownames(z_test[[i]]) = Entrez_exprs_test[[i]]$EntrezGene.ID
  colnames(z_test[[i]]) = colnames(df)
  z_test[[i]] = as.data.frame(z_test[[i]])
  z_test[[i]]$EntrezGene.ID = Entrez_exprs_test[[i]]$EntrezGene.ID
  rm(t, z_t, df)
}

names(z_test) = names(z)

# Joining in one expression matrix
# Test
z_test_exprs = z_test[[1]] %>% inner_join(z_test[[2]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[3]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[4]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[5]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[6]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[7]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[8]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[9]], by = "EntrezGene.ID") %>%
  inner_join(z_test[[10]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

z_test_exprs = na.omit(z_test_exprs)
ids = z_test_exprs$EntrezGene.ID
z_test_exprs = as.matrix(z_test_exprs %>% dplyr::select(-EntrezGene.ID))
z_test_exprs = z_test_exprs[, Pheno_test$Sample.ID]
rownames(z_test_exprs) = ids
z_test_exprs = z_test_exprs[rowSums(is.na(z_test_exprs)) != ncol(z_test_exprs), ]

# Training and Validation
z_exprs = z[[1]] %>% inner_join(z[[2]], by = "EntrezGene.ID") %>%
  inner_join(z[[3]], by = "EntrezGene.ID") %>%
  inner_join(z[[4]], by = "EntrezGene.ID") %>%
  inner_join(z[[5]], by = "EntrezGene.ID") %>%
  inner_join(z[[6]], by = "EntrezGene.ID") %>%
  inner_join(z[[7]], by = "EntrezGene.ID") %>%
  inner_join(z[[8]], by = "EntrezGene.ID") %>%
  inner_join(z[[9]], by = "EntrezGene.ID") %>%
  inner_join(z[[10]], by = "EntrezGene.ID") %>%
  dplyr::select(EntrezGene.ID, everything())

z_exprs = na.omit(z_exprs)
z_exprs = as.matrix(z_exprs %>% dplyr::select(-EntrezGene.ID))
z_exprs = z_exprs[, Pheno_exprs$Sample.ID]
rownames(z_exprs) = rownames(check_exprs)
z_exprs = z_exprs[rowSums(is.na(z_exprs)) != ncol(z_exprs), ]
KBZ_mds = plotMDS(z_exprs)
KBZ_pca = data.frame(cbind(KBZ_mds$x, KBZ_mds$y, as.character(Pheno_exprs$Dataset), Pheno_exprs$Sample.ID, 
                           Pheno_exprs$Treatment, as.character(Pheno_exprs$Response), 
                           as.character(Pheno_exprs$Timepoint_coded),
                           Pheno_exprs$scmod1, as.character(Pheno_exprs$pam50), 
                           as.character(Pheno_exprs$IC10), 
                           as.character(Pheno_exprs$Mammaprint_risk), 
                           as.character(Pheno_exprs$rorS_risk)))
KBZ_pca$Mammaprint = KBZ_pca$X11
KBZ_pca$Mammaprint[which(KBZ_pca$X11=="1")] = "Risk"
KBZ_pca$Mammaprint[which(KBZ_pca$X11!="1")] = "No risk"
KBZ_pca = KBZ_pca %>% dplyr::select(-X11)
colnames(KBZ_pca) = c("X1", "X2", "Dataset", "Sample.ID", "Treatment", "Response", "Timepoint",
                      "scmod1", "pam50", "IC10", "rorS", "Mammaprint")

for (i in c(3,5:12)){
  KBZ_pca[,i] = factor(KBZ_pca[,i], levels = unique(KBZ_pca[,i]), labels = unique(KBZ_pca[,i]))
}

class(KBZ_pca$X1) = "numeric"
class(KBZ_pca$X2) = "numeric"

z_pca = ggplot(KBZ_pca, aes(X1, X2, color = Dataset, shape = Response)) +
  geom_point(size = 0.2) +
  scale_color_brewer(palette = "Set3") +
  theme(plot.title = element_text(face = "bold", size = 5, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 size = 2.5),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 3.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.2),
        legend.position = "right",
        legend.key.size = unit(1, units = "mm"),
        legend.key.height = unit(1.5, "mm"),
        legend.text = element_text(size = 2.5),
        legend.title = element_text(face = "bold", size = 3),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(1, units = "mm"),
        legend.spacing.x = unit(1, units = "mm"))+
  labs(title = "Normalised samples - Multidimensional Scaling Plot",
       # x = paste0("\nMDS1 (", round(100*KBZ_mds$var.explained[1],2), "% of variance)"),
       # y = paste0("MDS2 (", round(100*KBZ_mds$var.explained[2],2), "% of variance)\n")
       x = "MDS1", y = "MDS2")
z_pca
ggsave(filename = "z_MDS.tiff",
       path = "QC", 
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Expression boxplots - Not helpful due to large number of samples #####
eset = exprs[,2:ncol(exprs)]
png("QC/Boxplot.png", width = 1920, height = 1080)
ggplot(melt(eset), aes(x=variable, y=value)) +
  geom_boxplot(outlier.size = 0.4, outlier.shape = 20,
               fill = c(rep("cyan", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C1",])), 
                        rep("orange", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C2",])),
                        rep("red4", nrow(Pheno_exprs[Pheno_exprs$Dataset=="C3",])),
                        rep("darkseagreen", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_1",])),
                        rep("beige", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_2",])),
                        rep("blue4", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E1_3",])),
                        rep("deeppink", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E2",])),
                        rep("purple", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E3",])),
                        rep("grey", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E4_1",])),
                        rep("darkgreen", nrow(Pheno_exprs[Pheno_exprs$Dataset=="E4_2",]))), 
               outlier.alpha = 0.1)+
  scale_y_continuous("Expression", limits = c(0,round(max(melt(eset)$value)+1)), 
                     breaks = seq(0,round(max(melt(eset)$value)+1), 1))+
  theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1,
                                   size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  size = 25, face = "bold"),
        axis.line = element_line())+
  labs(title = "Boxplot of expression",
       x = "\nSamples",
       y = "Expression\n")
dev.off()

# Heatmaps #####
save_pheatmap_png <- function(x, filename, width = 174, height = 120.8333, res = 650,
                              units = "mm") {
  png(filename, width = width, height = height, res = res, units = units, 
      compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_tiff <- function(x, filename, width = 174, height = 120.8333, res = 650,
                               units = "mm") {
  tiff(filename, width = width, height = height, res = res, units = units, 
       compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


annotation_for_heatmap = pca[, c(3, 5:12)]
rownames(annotation_for_heatmap) = pca$Sample.ID

dists = as.matrix(dist(t(plot_matrix), method = "manhattan"))

rownames(dists) = pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Response = c(Responder = "dodgerblue4", Non_responder = "deeppink4"),
  Dataset = c(C1 = "cyan", C2 = "orange", C3 = "red4", E1_1 = "darkseagreen",
              E1_2 = "beige", E1_3 = "blue4", E2 = "deeppink",
              E3 = "purple", E4_1 = "grey", E4_2 = "darkgreen"),
  Timepoint = c(T1 = "goldenrod2", T2 = "purple4"),
  Treatment = c(Chemotherapy = "paleturquoise3", Endocrine_treatment = "sienna2"),
  scmod1 = c(`ER-/HER2-` = "red4", `ER+/HER2- High Prolif` = "skyblue",
             `ER+/HER2- Low Prolif` = "darkblue", `HER2+` = "violet"),
  pam50 = c(Basal = "red4", LumB = "skyblue",
            LumA = "darkblue", Her2 = "purple4", Normal = "lightgreen"),
  IC10 = c(iC1 = "pink3", iC2 = "red", iC3 = "red4", iC4 = "beige", iC5 = "aliceblue",
           iC6 = "cadetblue3", iC7 = "darkmagenta", iC8 = "hotpink4", iC9 = "lightpink2",
           iC10 = "mistyrose1"),
  rorS = c(High = "red4", Intermediate = "orange", Low = "chartreuse4"),
  Mammaprint = c(Risk = "red3", `No risk` = "green")
)

heatmap = pheatmap(t(dists), col = (hmcol), 
                   annotation_col = annotation_for_heatmap,
                   annotation_colors = ann_colors,
                   cluster_cols = T,
                   treeheight_col =  0,
                   legend = TRUE,
                   show_rownames = F,
                   show_colnames = F,
                   fontsize = 3.5,
                   legend_breaks = c(min(dists, na.rm = TRUE), 
                                     max(dists, na.rm = TRUE)), 
                   legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(heatmap, "QC/Original_heatmap.tiff")

# RBE
RBE_annotation_for_heatmap = RBE_pca[, c(3, 5:12)]
rownames(RBE_annotation_for_heatmap) = RBE_pca$Sample.ID

rbe_dists = as.matrix(dist(t(RBE_exprs), method = "manhattan"))

rownames(rbe_dists) = RBE_pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(rbe_dists) <- NULL
diag(rbe_dists) <- NA

RBE_heatmap = pheatmap(t(rbe_dists), col = (hmcol), 
                       annotation_col = RBE_annotation_for_heatmap,
                       annotation_colors = ann_colors,
                       cluster_cols = T,
                       treeheight_col =  0,
                       legend = TRUE,
                       show_colnames = F,
                       show_rownames = F,
                       fontsize = 3.5,
                       legend_breaks = c(min(rbe_dists, na.rm = TRUE), 
                                         max(rbe_dists, na.rm = TRUE)), 
                       legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(RBE_heatmap, "QC/RBE_heatmap.tiff")

# KBZ
KBZ_annotation_for_heatmap = KBZ_pca[, c(3, 5:12)]
rownames(KBZ_annotation_for_heatmap) = KBZ_pca$Sample.ID

z_dists = as.matrix(dist(t(z_exprs), method = "manhattan"))

rownames(z_dists) = KBZ_pca$Sample.ID
hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
colnames(z_dists) <- NULL
diag(z_dists) <- NA

KBZ_heatmap = pheatmap(t(z_dists), col = (hmcol), 
                       annotation_col = KBZ_annotation_for_heatmap,
                       annotation_colors = ann_colors,
                       cluster_cols = T,
                       treeheight_col =  0,
                       legend = TRUE,
                       show_colnames = F,
                       show_rownames = F,
                       fontsize = 3.5,
                       legend_breaks = c(min(z_dists, na.rm = TRUE), 
                                         max(z_dists, na.rm = TRUE)), 
                       legend_labels = (c("small distance", "large distance")))
save_pheatmap_tiff(KBZ_heatmap, "QC/z_heatmap.tiff")

# 3D plots
library(rgl)
library(plot3D)
library(plot3Drgl)

exprs_for_join = exprs
str_sub(exprs_for_join$EntrezGene.ID, 0, 0) = "G_"
rownames(exprs_for_join) = exprs_for_join$EntrezGene.ID
exprs_for_join = as.data.frame(t(exprs_for_join))
exprs_for_join = exprs_for_join[2:nrow(exprs_for_join),]
exprs_for_join$Sample.ID = rownames(exprs_for_join)
masternormal = pca %>% inner_join(exprs_for_join, by =  "Sample.ID")

normaldist <- dist(masternormal[,13:ncol(masternormal)])
normalmds <- cmdscale(normaldist, k = 3)
PC1 <- normalmds[,1]
PC2 <- normalmds[,2]
PC3 <- normalmds[,3]

#Plot
scatter3D(PC1, PC2, PC3, phi = 10, bty = "u", theta = 45, width = 1920, height = 1080,
          main = "3D PCA - raw data", xlab = "PC1", ylab = "PC2", zlab = "PC3",
          ticktype = "detailed", pch = 20, alpha = 0.7,
          colvar = as.integer(masternormal$Dataset),
          col = c("cyan", "orange", "red4", "darkseagreen",
                  "beige", "blue4", "deeppink",
                  "purple", "grey", "darkgreen"), 
          colkey = list(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                        addlines = TRUE, length = 1, width = 1,
                        labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3",
                                   "E2", "E3", "E4_1", "E4_2")), clab = "Dataset", cex = 1)
plotrgl()

# KBZ #
z_exprs_for_join = as.data.frame(t(z_exprs))
str_sub(colnames(z_exprs_for_join), 0, 0) = "G_"
z_exprs_for_join$Sample.ID = rownames(z_exprs_for_join)
KBZ_master = KBZ_pca %>% inner_join(z_exprs_for_join, by =  "Sample.ID")

KBZdist <- dist(KBZ_master[,13:ncol(KBZ_master)])
KBZmds <- cmdscale(KBZdist, k = 3)
KBZPC1 <- KBZmds[,1]
KBZPC2 <- KBZmds[,2]
KBZPC3 <- KBZmds[,3]

# Plot
scatter3D(KBZPC1, KBZPC2, KBZPC3, phi = 10, bty = "u", theta = 45, width = 1920, height = 1080,
          main = "3D PCA - z-scores", xlab = "PC1", ylab = "PC2", zlab = "PC3",
          ticktype = "detailed", pch = 20, alpha = 0.7,
          colvar = as.integer(KBZ_master$Dataset),
          col = c("cyan", "orange", "red4", "darkseagreen",
                  "beige", "blue4", "deeppink",
                  "purple", "grey", "darkgreen"), 
          colkey = list(at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                        addlines = TRUE, length = 1, width = 1,
                        labels = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3",
                                   "E2", "E3", "E4_1", "E4_2")), clab = "Dataset", cex = 1)
plotrgl()

rm(i, j, k)

##### Plotly #####
library(plotly)
library(data.table)
library(colorspace)

Pheno_sunburst = Pheno_full %>% dplyr::select(Treatment, Dataset, Response, pam50)

Pheno_sunburst = Pheno_sunburst %>%
  group_by(Treatment, Dataset, Response, pam50) %>%
  summarise(Counts = n()) %>%
  as.data.frame()

as.sunburstDF <- function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF <- names(DF)
  
  if(is.data.table(DF)){
    DT <- copy(DF)
  } else {
    DT <- data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT <- names(DT)
  hierarchy_columns <- setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList <- list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns <- colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT <- unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT <- DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <- rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns <- setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}

coloring  = data.frame(stringsAsFactors = FALSE,
                       colors = tolower(gplots::col2hex(c("paleturquoise3", "sienna2", "paleturquoise1",
                                                          "paleturquoise1", "paleturquoise1", "sienna1", 
                                                          "sienna1", "sienna1", "sienna1",
                                                          "sienna1", "sienna1", "sienna1", "deeppink4", "dodgerblue4", 
                                                          "red4", "violet", "darkblue", "skyblue", "lightgreen"))),
                       labels = c("Chemotherapy", "Endocrine_treatment", "C1", "C2", "C3", 
                                  "E1_1", "E1_2", "E1_3", "E2", "E3", "E4_1", "E4_2",
                                  "Non_responder", "Responder", "Basal", "Her2", "LumA", 
                                  "LumB", "Normal"))

sunburstDF <- as.sunburstDF(Pheno_sunburst, value_column = "Counts", add_root = FALSE) %>%
  inner_join(coloring, by = "labels")

pie = plot_ly() %>%
  add_trace(ids = sunburstDF$ids, labels= sunburstDF$labels, parents = sunburstDF$parents, 
            values= sunburstDF$values, type='sunburst', branchvalues = 'total',
            insidetextorientation='radial', maxdepth = 4,
            marker = list(colors = sunburstDF$colors)) %>%
  layout(
    grid = list(columns =1, rows = 1),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )
pie

# Session Info #####
sessionInfo()
>>>>>>> 494bac1638efd88b01d931820c8951ad15a7290e
