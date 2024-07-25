### This script can be used to make files for each participant and visit to upload in PEP, 
### with the results of blood inflammation markers (CRP, NFL, Olink and COVID antibodies)

rm(list = ls())

#### Load all raw data files, one for CRP/NFL, one for Olink and one for COVID antibodies ####

library("readxl")
library(dplyr)

# CRP and NFL file
CRPNFL <- read_excel("M:/Documents/Projecten/Inflammatie/Raw_data_files/Raw_Result_PPP_NFL_CRP.xlsx") %>%
                     select(Sampl, ItemID, Datum, NFL, CRP, Visit) %>%
                     filter(!Sampl %in% c("POM1PL3107813", "POM2PL0755946", "POM3PL4994453")) %>%  # No sample after COVID, not included in measurements so remove
                     mutate(Datum  = as.Date(Datum),
                            NFL    = as.double(gsub("[^0-9.-]", "", NFL)), 
                            CRP    = as.double(gsub("\\*", "", gsub("<0.17", "0.10", CRP))),       # values <LOD are not excluded but counted as 0.10 (between 0-0.17), * before values are removed
                            Sampl  = as.factor(Sampl),
                            ItemID = as.factor(ItemID),
                            Visit  = as.factor(Visit))
CRPNFL$Datum <- format(CRPNFL$Datum, "%d-%h-%Y")

# Covid antibodies; column names in original file are long and with spaces, so in sheet 2 I adjusted the names for loading in R
Cov_anti <- read_excel("M:/Documents/Projecten/Inflammatie/Raw_data_files/Raw_Result_PPP_Covid_antibodies.xlsx", sheet = 2) %>%
                       mutate(ItemID       = as.factor(ItemID),
                              N_protein    = as.double(N_protein), # NAs are introduced, should later be replaced by OOR<
                              S_protein    = as.double(S_protein), 
                              RBD_WuhanHu = as.double(RBD_WuhanHu),
                              RBD_delta    = as.double(RBD_delta),
                              RBD_omicron  = as.double(RBD_omicron))                        

# Olink marker file; use sheet 2 where the names are adjusted for loading in R
Olink <- read_excel("M:/Documents/Projecten/Inflammatie/Raw_data_files/Raw_Result_PPP_Olink_NPX.xlsx", sheet = 2) %>%
                    mutate(ItemID = as.factor(ItemID))

# Make sure that all files contain sample ID, date and visit number; use ItemID to merge (this variable is present in all files)
Cov_anti <- merge(Cov_anti, CRPNFL[, c("ItemID", "Sampl", "Visit", "Datum")], by = "ItemID", all.x = TRUE) 
Olink <- merge(Olink, CRPNFL[, c("ItemID", "Sampl", "Visit", "Datum")], by = "ItemID", all.x = TRUE) 

#### Make files for every sample to upload in PEP for CRP and NfL ####

dirUpload <- "M:/Documents/Projecten/Inflammatie/PEP_upload_files"

subfolder_names <- c("NFL.Visit1", "NFL.Visit2", "NFL.Visit3", 
                     "CRP.Visit1", "CRP.Visit2", "CRP.Visit3",
                     "NFL.Visit1.Date", "NFL.Visit2.Date", "NFL.Visit3.Date", 
                     "CRP.Visit1.Date", "CRP.Visit2.Date", "CRP.Visit3.Date")

for (s in subfolder_names) {
  dir.create(file.path(dirUpload, s), recursive = TRUE)
}

for (i in 1:nrow(CRPNFL)) {
  ID <- CRPNFL$Sampl[i]
  
  if (CRPNFL$Visit[i] == "Visit 1") {
    
    # Concentration files
    cFile_NFL <- file.path(dirUpload, "NFL.Visit1", paste0(ID, "_InflammationMarkers_Blood_NFL.Visit1.tsv"))
    cNFL <- CRPNFL$NFL[i] # Make one row for NFL
    write.table(cNFL, cFile_NFL, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    cFile_CRP <- file.path(dirUpload, "CRP.Visit1", paste0(ID, "_InflammationMarkers_Blood_CRP.Visit1.tsv"))
    cCRP <- CRPNFL$CRP[i]
    write.table(cCRP, cFile_CRP, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date files
    cFile_NFL_date <- file.path(dirUpload, "NFL.Visit1.Date", paste0(ID, "_InflammationMarkers_Blood_Date_NFL.Visit1.tsv"))
    cFile_CRP_date <- file.path(dirUpload, "CRP.Visit1.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CRP.Visit1.tsv"))
    cDate <- CRPNFL$Datum[i]
    write.table(cDate, cFile_NFL_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(cDate, cFile_CRP_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (CRPNFL$Visit[i] == "Visit 2") {
    
    # Concentration files
    cFile_NFL <- file.path(dirUpload, "NFL.Visit2", paste0(ID, "_InflammationMarkers_Blood_NFL.Visit2.tsv"))
    cNFL <- CRPNFL$NFL[i] # Make one row for NFL
    write.table(cNFL, cFile_NFL, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    cFile_CRP <- file.path(dirUpload, "CRP.Visit2", paste0(ID, "_InflammationMarkers_Blood_CRP.Visit2.tsv"))
    cCRP <- CRPNFL$CRP[i]
    write.table(cCRP, cFile_CRP, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date files
    cFile_NFL_date <- file.path(dirUpload, "NFL.Visit2.Date", paste0(ID, "_InflammationMarkers_Blood_Date_NFL.Visit2.tsv"))
    cFile_CRP_date <- file.path(dirUpload, "CRP.Visit2.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CRP.Visit2.tsv"))
    cDate <- CRPNFL$Datum[i]
    write.table(cDate, cFile_NFL_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(cDate, cFile_CRP_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (CRPNFL$Visit[i] == "Visit 3") {
    
    # Concentration files
    cFile_NFL <- file.path(dirUpload, "NFL.Visit3", paste0(ID, "_InflammationMarkers_Blood_NFL.Visit3.tsv"))
    cNFL <- CRPNFL$NFL[i]  # Make one row for NFL
    write.table(cNFL, cFile_NFL, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    cFile_CRP <- file.path(dirUpload, "CRP.Visit3", paste0(ID, "_InflammationMarkers_Blood_CRP.Visit3.tsv"))
    cCRP <- CRPNFL$CRP[i]
    write.table(cCRP, cFile_CRP, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date files
    cFile_NFL_date <- file.path(dirUpload, "NFL.Visit3.Date", paste0(ID, "_InflammationMarkers_Blood_Date_NFL.Visit3.tsv"))
    cFile_CRP_date <- file.path(dirUpload, "CRP.Visit3.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CRP.Visit3.tsv"))
    cDate <- CRPNFL$Datum[i]
    write.table(cDate, cFile_NFL_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(cDate, cFile_CRP_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  }
}

#### Make files for every sample to upload in PEP for Covid antibodies #### 

subfolder_names2 <- c("CovidAntibodies.Visit2", "CovidAntibodies.Visit3", "CovidAntibodies.Visit2.Date", "CovidAntibodies.Visit3.Date") # There are no samples from visit 1 for COVID antibodies

for (s in subfolder_names2) {
  dir.create(file.path(dirUpload, s), recursive = TRUE)
}

for (i in 1:nrow(Cov_anti)) {
  ID <- Cov_anti$Sampl[i]
  
  if (Cov_anti$Visit[i] == "Visit 1") {
    
    # Concentration file
    cFile_Covid <- file.path(dirUpload, "CovidAntibodies.Visit1", paste0(ID, "_InflammationMarkers_Blood_CovidAntibodies.Visit1.tsv"))
    cCovid <- data.frame(Name = c("S protein", "N protein", "RBD Wuhan Hu 1", "RBD delta", "RBD omicron"),
                         Value = c(Cov_anti$S_protein[i], Cov_anti$N_protein[i], Cov_anti$RBD_WuhanHu[i], Cov_anti$RBD_delta[i], Cov_anti$RBD_omicron[i]))
    write.table(cCovid, cFile_Covid, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "CovidAntibodies.Visit1.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CovidAntibodies.Visit1.tsv"))
    cDate <- Cov_anti$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (Cov_anti$Visit[i] == "Visit 2") {
    
    # Concentration file
    cFile_Covid <- file.path(dirUpload, "CovidAntibodies.Visit2", paste0(ID, "_InflammationMarkers_Blood_CovidAntibodies.Visit2.tsv"))
    cCovid <- data.frame(Name = c("S protein", "N protein", "RBD Wuhan Hu 1", "RBD delta", "RBD omicron"),
                         Value = c(Cov_anti$S_protein[i], Cov_anti$N_protein[i], Cov_anti$RBD_WuhanHu[i], Cov_anti$RBD_delta[i], Cov_anti$RBD_omicron[i]))
    write.table(cCovid, cFile_Covid, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "CovidAntibodies.Visit2.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CovidAntibodies.Visit2.tsv"))
    cDate <- Cov_anti$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (Cov_anti$Visit[i] == "Visit 3") {
    
    # Concentration file
    cFile_Covid <- file.path(dirUpload, "CovidAntibodies.Visit3", paste0(ID, "_InflammationMarkers_Blood_CovidAntibodies.Visit3.tsv"))
    cCovid <- data.frame(Name = c("S protein", "N protein", "RBD Wuhan Hu 1", "RBD delta", "RBD omicron"),
                         Value = c(Cov_anti$S_protein[i], Cov_anti$N_protein[i], Cov_anti$RBD_WuhanHu[i], Cov_anti$RBD_delta[i], Cov_anti$RBD_omicron[i]))
    write.table(cCovid, cFile_Covid, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "CovidAntibodies.Visit3.Date", paste0(ID, "_InflammationMarkers_Blood_Date_CovidAntibodies.Visit3.tsv"))
    cDate <- Cov_anti$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  }
}


#### Make files for every sample to upload in PEP for Olink #### 

subfolder_names3 <- c("Olink96.Visit1", "Olink96.Visit2", "Olink96.Visit3", 
                      "Olink96.Visit1.Date", "Olink96.Visit2.Date", "Olink96.Visit3.Date")

all_markers <- 2:97 # Olink marker concentrations and QC info are in columns 2:97 in the adjusted sheet of the xlsx file

for (s in subfolder_names3) {
  dir.create(file.path(dirUpload, s), recursive = TRUE)
}

for (i in 1:nrow(Olink)) {
  ID <- Olink$Sampl[i]
  
  if (Olink$Visit[i] == "Visit 1") {
    
    # Concentration file
    cFile_Olink <- file.path(dirUpload, "Olink96.Visit1", paste0(ID, "_InflammationMarkers_Blood_Olink96.Visit1.tsv"))
    cOlink <- data.frame(Olink[i, all_markers])
    write.table(cOlink, cFile_Olink, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "Olink96.Visit1.Date", paste0(ID, "_InflammationMarkers_Blood_Date_Olink96.Visit1.tsv"))
    cDate <- Olink$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (Olink$Visit[i] == "Visit 2") {
    
    # Concentration file
    cFile_Olink <- file.path(dirUpload, "Olink96.Visit2", paste0(ID, "_InflammationMarkers_Blood_Olink96.Visit2.tsv"))
    cOlink <- data.frame(Olink[i, all_markers])
    write.table(cOlink, cFile_Olink, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "Olink96.Visit2.Date", paste0(ID, "_InflammationMarkers_Blood_Date_Olink96.Visit2.tsv"))
    cDate <- Olink$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)
    
  } else if (Olink$Visit[i] == "Visit 3") {
    
    # Concentration file
    cFile_Olink <- file.path(dirUpload, "Olink96.Visit3", paste0(ID, "_InflammationMarkers_Blood_Olink96.Visit3.tsv"))
    cOlink <- data.frame(Olink[i, all_markers])
    write.table(cOlink, cFile_Olink, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    # Date file
    cFile_date <- file.path(dirUpload, "Olink96.Visit3.Date", paste0(ID, "_InflammationMarkers_Blood_Date_Olink96.Visit3.tsv"))
    cDate <- Olink$Datum[i]
    write.table(cDate, cFile_date, sep = "\t", row.names = FALSE, col.names = FALSE)    
  }
}
