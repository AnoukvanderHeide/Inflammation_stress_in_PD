# LET OP: in de huidige PEP download missen waarden boven reference value, dit betreft ~6 waarden van visit 1, die zijn bij de wide_adj file aangevuld

rm(list = ls())

# Load necessary library
library(data.table)
library(lubridate)
library(dplyr)
library(tidyr)

##### Load data of CRP and NFL concentrations in one dataframe #####

# Function to read data from the folders and create data frames per variable
read_data <- function(folder_path, measure) {
  folder_list <- list.files(path = folder_path, pattern = "^POMU", full.names = TRUE)
  data <- data.frame(Subj_ID = character(0), Value = numeric(0))
  
  for (folder in folder_list) {
    ID <- basename(folder)
    file_path <- file.path(folder, paste0("DD_InflammationMarkers_Blood_", measure))
    value <- as.numeric(readLines(file_path))
    data <- rbind(data, data.frame(Subj_ID = ID, Value = value))
  }
  colnames(data)[2] <- measure 
  return(data)
}

# Specify the folder paths and measures
folder_paths <- c(CRP.Visit1 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit1",
                  CRP.Visit2 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit2",
                  CRP.Visit3 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit3",
                  NFL.Visit1 = "M:/Documents/Projecten/Inflammatie/PEP_download/NFL.Visit1",
                  NFL.Visit2 = "M:/Documents/Projecten/Inflammatie/PEP_download/NFL.Visit2",
                  NFL.Visit3 = "M:/Documents/Projecten/Inflammatie/PEP_download/NFL.Visit3")

# Read and merge data for each measure
data_list <- lapply(names(folder_paths), function(measure) read_data(folder_paths[measure], measure))
Data_CRPNFL_wide <- Reduce(function(x, y) merge(x, y, by = "Subj_ID", all.x = TRUE), data_list)

# For the dates of the blood samples, different function is needed for _Date in folder name
read_data_dates <- function(folder_path, measure) {
  folder_list <- list.files(path = folder_path, pattern = "^POMU", full.names = TRUE)
  data <- data.frame(Subj_ID = character(0), Date_str = character(0))
  
  for (folder in folder_list) {
    ID <- basename(folder)
    file_path <- file.path(folder, paste0("DD_InflammationMarkers_Blood_Date_", measure))
    date_str <- readLines(file_path)
    data <- rbind(data, data.frame(Subj_ID = ID, Date_str = date_str))
  }
  colnames(data)[2] <- measure 
  return(data)
}

# Specify the folder paths and measures
folder_paths <- c(CRP.Visit1 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit1.Date",
                  CRP.Visit2 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit2.Date",
                  CRP.Visit3 = "M:/Documents/Projecten/Inflammatie/PEP_download/CRP.Visit3.Date")

# Read and merge data for each measure
data_list <- lapply(names(folder_paths), function(measure) read_data_dates(folder_paths[measure], measure))
Dates <- Reduce(function(x, y) merge(x, y, by = "Subj_ID", all.x = TRUE), data_list)

# Change to corresponding date format
date_columns <- colnames(Dates)[grepl("CRP.", colnames(Dates))]
Dates[, date_columns] <- lapply(Dates[, date_columns], dmy)     # Assumes lubridate's dmy function
colnames(Dates) <- gsub("CRP.", "Date.", colnames(Dates))       # Rename date columns 
Data_CRPNFL_wide <- merge(Data_CRPNFL_wide, Dates, by = "Subj_ID", all.x = TRUE) 
Data_CRPNFL_wide$Subj_ID <- paste("sub-", Data_CRPNFL_wide$Subj_ID, sep = "")

# Reshape the wide data frame to long format
Data_CRPNFL_long <- pivot_longer(Data_CRPNFL_wide, 
                                 cols = starts_with("CRP.Visit") | starts_with("NFL.Visit") | starts_with("Date.Visit"), 
                                 names_to = c(".value", "Visit"),
                                 names_pattern = "(CRP|NFL|Date)\\.Visit(\\d+)",
                                 values_to = c("CRP", "NFL", "Date"))
Data_CRPNFL_long <- Data_CRPNFL_long %>% filter(!is.na(CRP) | !is.na(NFL))

# sub-POMU1CF55A8FB405D10A	visit one and two are a duplicate, the values are from visit 2 so remove visit 1
Data_CRPNFL_wide <- Data_CRPNFL_wide %>%
  mutate(across(matches("Visit1"), ~ ifelse(Subj_ID == "sub-POMU1CF55A8FB405D10A", NA, .)))
Data_CRPNFL_long <- Data_CRPNFL_long %>%
  filter(!(Subj_ID == "sub-POMU1CF55A8FB405D10A" & Visit == "1"))

# Save files
write.csv(Data_CRPNFL_wide, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_CRPNFL_wide.csv", row.names = FALSE)
write.csv(Data_CRPNFL_long, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_CRPNFL_long.csv", row.names = FALSE)

##### For COVID antibodies #####

# New function for covid antibodies, where there are multiple markers per file
read_data_covid <- function(folder_path, measure) {
  folder_list <- list.files(path = folder_path, pattern = "^POMU", full.names = TRUE)
  data <- data.frame(Subj_ID = character(0), stringsAsFactors = FALSE)
  
  for (folder in folder_list) {
    ID <- basename(folder)
    file_path <- file.path(folder, paste0("DD_InflammationMarkers_Blood_", measure))
    lines <- readLines(file_path)
    values <- as.numeric(sub(".*\t", "", lines))
    col_names <- sub("\t.*", "", lines)
    row_data <- c(ID, values)
    
    # Create or extend the data frame
    data <- rbind(data, row_data)
    colnames(data) <- c("Subj_ID", col_names)
  }
  
  return(data)
}

# Read data for COVID antibodies for both visit 2 and 3
folder_visit2 <- "M:/Documents/Projecten/Inflammatie/PEP_download/CovidAntibodies.Visit2"
Covid_visit2 <- read_data_covid(folder_visit2, "CovidAntibodies.Visit2")
folder_visit3 <- "M:/Documents/Projecten/Inflammatie/PEP_download/CovidAntibodies.Visit3"
Covid_visit3 <- read_data_covid(folder_visit3, "CovidAntibodies.Visit3")

# Change column names and add visit column so that they are the same; this way they can be combined in full table.
colnames(Covid_visit2) <- c("Subj_ID", "S_protein", "N_protein", "RBD_Wuhan", "RBD_delta", "RBD_omicron")
colnames(Covid_visit3) <- c("Subj_ID", "S_protein", "N_protein", "RBD_Wuhan", "RBD_delta", "RBD_omicron")
Covid_visit2$Visit_Covid <- 2 # Specify at which visit the Covid antibodies were measured
Covid_visit3$Visit_Covid <- 3 # Specify at which visit the Covid antibodies were measured

# Convert marker values to numeric
marker_columns <- colnames(Covid_visit2)[-1]  # Exclude the first column (Subj_ID)
Covid_visit2[, marker_columns] <- sapply(Covid_visit2[, marker_columns], function(x) as.numeric(gsub(",", "", x)))
Covid_visit3[, marker_columns] <- sapply(Covid_visit3[, marker_columns], function(x) as.numeric(gsub(",", "", x)))

# Merge data for both visits
Data_Covid <- rbind(Covid_visit2, Covid_visit3)
Data_Covid$Subj_ID <- paste("sub-", Data_Covid$Subj_ID, sep = "")
Data_Covid$Visit_Covid <- as.factor(Data_Covid$Visit_Covid)

duplicated_rows <- Data_Covid[duplicated(Data_Covid$Subj_ID), ] # Check whether there are duplicates in subj ID (not the case)
Data_Covid <- Data_Covid[order(Data_Covid$Subj_ID), ]
Data_CRPNFL_wide <- Data_CRPNFL_wide[order(Data_CRPNFL_wide$Subj_ID), ]
Data_Covid <- Data_Covid %>% mutate(Date = case_when(Visit_Covid == 2 ~ Data_CRPNFL_wide$Date.Visit2,
                                                     Visit_Covid == 3 ~ Data_CRPNFL_wide$Date.Visit3,
                                                     TRUE ~ NA))

write.csv(Data_Covid, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_Covid_antibodies.csv", row.names = FALSE)


##### For Olink markers #####

# New function to extract Olink marker concentrations



read_data_olink96 <- function(folder_path, visit) {
  folder_list <- list.files(path = folder_path, pattern = "^POMU", full.names = TRUE)
  data <- data.frame(Subj_ID = character(0), stringsAsFactors = FALSE)
  
  for (folder in folder_list) {
    ID <- basename(folder)
    file_path <- file.path(folder, paste0("DD_InflammationMarkers_Blood_Olink96.", visit)) # Specify correct file
    lines <- readLines(file_path) # Load file
    col_names <- unlist(strsplit(lines[1], "\t"))
    values <- matrix(unlist(strsplit(lines[2], "\t")), ncol = length(col_names), byrow = TRUE)
    row_data <- c(ID, values)
    data <- rbind(data, row_data) # Create or extend the data frame
    colnames(data) <- c("Subj_ID", col_names) # Add marker names as column names 
  }
  
  return(data)
}


# Specify the folder paths and measures
folder_olink96_visit1 <- "M:/Documents/Projecten/Inflammatie/PEP_download/Olink96.Visit1"
Olink96_visit1 <- read_data_olink96(folder_olink96_visit1, "Visit1")
folder_olink96_visit2 <- "M:/Documents/Projecten/Inflammatie/PEP_download/Olink96.Visit2"
Olink96_visit2 <- read_data_olink96(folder_olink96_visit2, "Visit2")
folder_olink96_visit3 <- "M:/Documents/Projecten/Inflammatie/PEP_download/Olink96.Visit3"
Olink96_visit3 <- read_data_olink96(folder_olink96_visit3, "Visit3")

# Specify at which visit the Olink markers were assessed
Olink96_visit1$Visit <- 1 
Olink96_visit2$Visit <- 2 
Olink96_visit3$Visit <- 3

# Combine the data from all visits
Data_Olink96 <- rbind(Olink96_visit1, Olink96_visit2, Olink96_visit3)
colnames(Data_Olink96) <- gsub("\"", "", colnames(Data_Olink96))

# Add 'sub-' prefix to Subj_ID
Data_Olink96$Subj_ID <- paste("sub-", Data_Olink96$Subj_ID, sep = "")

# Add visit number to Subj_ID
Data_Olink96$Subj_ID <- ifelse(Data_Olink96$Visit == 1, paste(Data_Olink96$Subj_ID, "_vis1", sep = ""),
                               ifelse(Data_Olink96$Visit == 2, paste(Data_Olink96$Subj_ID, "_vis2", sep = ""),
                                      ifelse(Data_Olink96$Visit == 3, paste(Data_Olink96$Subj_ID, "_vis3", sep = ""),
                                             Data_Olink96$Subj_ID)))
Data_Olink96$Visit <- NULL # Remove visit column, since visit number is now included in subject ID
Data_Olink96 <- Data_Olink96[!(Data_Olink96$Subj_ID == "sub-POMU1CF55A8FB405D10A_vis1"),] # For some reason this one is duplicated in PEP and loaded for visit 1 and visit 2 both


write.csv(Data_Olink96, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink96_Routput.csv", row.names = FALSE)
