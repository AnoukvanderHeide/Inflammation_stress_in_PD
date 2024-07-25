
rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)

##### -------- COVID antibodies concentrations: identify vaccinated and infected people from blood samples -------- #####

Data_Covid <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_Covid_antibodies.csv")

Data_Covid$S_protein_status <- ifelse(Data_Covid$S_protein > 2.28, "+", "-") # Use average for S-protein in neg. controls as cutoff
Data_Covid$N_protein_status <- ifelse(Data_Covid$N_protein > 8.21, "+", "-") # Use average for N-protein in neg. controls as cutoff
Data_Covid$N_protein_log <- log(Data_Covid$N_protein)
Data_Covid$S_protein_log <- log(Data_Covid$S_protein)

print(xtabs(~ S_protein_status + N_protein_status, data = Data_Covid))

# Determine COVID infection status at time of blood sample
Data_Covid$infected <- as.factor(ifelse(Data_Covid$Date < as.Date("2020-03-01"), "precov", # there is no sample pre-COVID but just to check
                                        ifelse(Data_Covid$N_protein_status == "+", "yes", "no")))
table(Data_Covid$infected)

# Plot all values N-protein (log-transformed)
ggplot(Data_Covid, aes(x = 1, y = N_protein_log)) +
        geom_jitter(aes(color = infected), alpha = 0.6, na.rm = TRUE) +
        geom_hline(yintercept = log(8.21), linetype = "dashed", color = "red") +  # Add horizontal line
        labs(x = "", y = "ln(N-protein)") + 
        theme_minimal() +
        theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))

# Determine vaccination status at time of blood sample
Data_Covid$prevac <- ifelse(Data_Covid$Date < as.Date("2021-03-01"), 1, 0) # 184 were before vaccination 
Data_Covid$vaccinated <- as.factor(ifelse(Data_Covid$Date < as.Date("2021-03-01"), "prevac", 
                                          ifelse(Data_Covid$S_protein_status == "+" & Data_Covid$N_protein_status == "-", "yes", 
                                                ifelse(Data_Covid$S_protein_status == "+" & Data_Covid$N_protein_status == "+", "possible", "probably unvaccinated"))))
filtered_postvac <- Data_Covid %>% filter(prevac == 0)
table(Data_Covid$vaccinated)

# all values S-protein (log-transformed)
ggplot(filtered_postvac, aes(x = 1, y = S_protein_log)) +
        geom_jitter(aes(color = vaccinated), alpha = 0.6) +  # Adjust width and alpha for jitter
        geom_hline(yintercept = log(2.28), linetype = "dashed", color = "red") +  # Add horizontal line
        labs(x = "", y = "ln(S-protein)") +
        theme_minimal() +
        theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))

# only values <50
#filtered_low <- filtered_postvac %>% filter(S_protein < 50)
#ggplot(filtered_low, aes(x = 1, y = S_protein)) +
#  geom_jitter(alpha = 0.6) +  # Adjust width and alpha for jitter
#  geom_hline(yintercept = 2.28, linetype = "dashed", color = "red") +  # Add horizontal line
#  labs(x = "", y = "S_protein Values")

# separate panels S-protein
#ggplot(filtered_postvac, aes(x = 1, y = S_protein)) +
#  geom_jitter(aes(color = vaccinated), width = 0.1, height = 0, alpha = 0.6) +
#  geom_hline(yintercept = 2.28, linetype = "dashed") +
#  facet_wrap(~ifelse(!is.na(S_protein) & S_protein < 50, "Below 50", "Above 50"), scales = "free_y", ncol = 1) +
#  labs(x = "", y = "S_protein concentration")

write.csv(Data_Covid, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_Covid_analysed.csv", row.names = FALSE)


##### -------- COVID antibodies: identify infected people from survey responses -------- #####

subj_numbers <- Data_Covid$Subj_ID
directory_path <- "P:/3022026.01/pep/ClinVars4"
library(jsonlite)

### COVID biweekly survey responses ###

# Initialize an empty list to store df for each subject
Biweekly_list <- list()

# Loop over all subjects
for (subj in subj_numbers) {
  subj_directory  <- file.path(directory_path, paste0(subj))
  
  if (file.exists(subj_directory)) {
    json_files  <- list.files(path = file.path(subj_directory, "ses-COVIDweek2"),
                              pattern = "COVID.Castor.CovidQuestionnaires.CovPackWeek2.CovDiagno.CovDiagno.AnswerSet.*.json",
                              full.names = TRUE)
    json_files  <- json_files[!grepl("WeekNumber", json_files)]    # Filter files without WeekNumber in their names
    subj_values <- list()                                          # Make a list to store values for this subject
    
    # Loop through file paths
    for (i in 0:8) {
      file_path <- json_files[grepl(paste0(".AnswerSet", i, ".json"), json_files)]
      if (length(file_path) == 0) {subj_values[[paste0("AnswerSet", i)]] <- NA} 
      else {json_data    <- fromJSON(file_path)
      contaminated <- json_data$crf$ContaminatedLstWeek1    # is 1 if Cov19TstResLstWeek = 1 OR Cov19DoctAssessLstWeek = 1 OR Cov19SelfAssessLstWeek = 1
      subj_values[[paste0("AnswerSet", i)]] <- contaminated # Store values in the list
      }
    }
    
    Biweekly_subject <- data.frame(Subj_ID = subj, do.call(cbind, subj_values)) # Create a data frame for this subject
    Biweekly_list[[length(Biweekly_list) + 1]] <- Biweekly_subject} 
  else {cat(sprintf("Subject directory for %s does not exist.\n", subj))
  }
}

# Combine data from all subjects
Biweekly <- do.call(rbind, Biweekly_list)
Biweekly[, 2:ncol(Biweekly)] <- sapply(Biweekly[, 2:ncol(Biweekly)], as.numeric)

### COVID baseline and final survey responses ###

# Create an empty dataframe
Baseline <- data.frame(Subj_ID = subj_numbers,
                       Contam1 = NA, Contam2 = NA, Selftest = NA, Covidtest = NA, Doctor = NA)

# Loop over subjects and store responses
for (subj in subj_numbers) {
  subj_directory        <- file.path(directory_path, paste0(subj))
  if (file.exists(subj_directory)) {
    file_path         <- file.path(subj_directory, "ses-COVIDbasic", "COVID.Castor.CovidQuestionnaires.CovPackBasic.CovComDia.CovComDia.AnswerSet0.json")
    
    if (file.exists(file_path)) {
      json_data     <- fromJSON(file_path)
      contaminated1 <- json_data$crf$Contaminated1 # Cov19TstRes = 1 OR Cov19DoctAssess = 1 OR Cov19SelfAssess = 1
      contaminated2 <- json_data$crf$Contaminated2
      selftest      <- json_data$crf$Cov19SelfAssess
      covidtest     <- json_data$crf$CovidTest
      doctor        <- json_data$crf$Cov19DoctAssess
      
      # Find the row index where Subj_ID matches and update the values
      row_index <- which(Baseline$Subj_ID == subj)
      Baseline[row_index, "Contam1"] <- as.numeric(contaminated1)
      Baseline[row_index, "Contam2"] <- as.numeric(contaminated2)
      Baseline[row_index, "Selftest"] <- as.numeric(selftest)
      Baseline[row_index, "Covidtest"] <- as.numeric(covidtest)
      Baseline[row_index, "Doctor"] <- as.numeric(doctor)} 
    else {cat(sprintf("JSON file for subject %s does not exist.\n", subj))
    }} 
  else {cat(sprintf("Subject directory for %s does not exist.\n", subj))
  }
}

sum(Baseline$Contam1 == 1, na.rm = TRUE)   # 17 people contaminated at baseline survey (not sure what that means)
sum(Baseline$Selftest == 1, na.rm = TRUE)  # 17 people infection determined by selftest (of self-assessment? tests not yet available?)
sum(Baseline$Covidtest == 1, na.rm = TRUE) # for 1 infection determined by covid-test
sum(Baseline$Doctor == 1, na.rm = TRUE)    # for 2 infection determined by doctor


# Create an empty dataframe for final survey responses
Final <- data.frame(Subj_ID = subj_numbers, Contam11 = NA)

# Loop over subjects and store responses
for (subj in subj_numbers) {
  subj_directory <- file.path(directory_path, paste0(subj))
  if (file.exists(subj_directory)) {
    file_path <- file.path(subj_directory, "ses-COVIDfinal", "COVID.Castor.CovidQuestionnaires.CovPackFinal.CovDiagno.CovDiagno.AnswerSet0.json")
    if (file.exists(file_path)) {
      json_data <- fromJSON(file_path)
      Contam11 <- json_data$crf$ContaminatedLstWeek1 # Cov19TstRes = 1 OR Cov19DoctAssess = 1 OR Cov19SelfAssess = 1
      row_index <- which(Final$Subj_ID == subj)      # Find the row index where Subj_ID matches and update the values
      Final[row_index, "Contam11"] <- as.numeric(Contam11)} 
    else {cat(sprintf("JSON file for subject %s does not exist.\n", subj))
    }}
  else {cat(sprintf("Subject directory for %s does not exist.\n", subj))
  }
}

### COVID combined survey responses ###

All_surveys <- merge(Biweekly, Baseline[, c("Subj_ID", "Contam1")], by = "Subj_ID", all.x = TRUE)
All_surveys <- merge(All_surveys, Final[, c("Subj_ID", "Contam11")], by = "Subj_ID", all.x = TRUE)
colnames(All_surveys)[colnames(All_surveys) %in% paste0("AnswerSet", 0:8)] <- paste0("Contam", 2:10)
All_surveys <- All_surveys %>% relocate("Contam1", .before="Contam2")
All_surveys <- All_surveys %>% relocate("Contam11", .after="Contam10")

# Check nr. of infections according to surveys
library(dplyr)
All_surveys <- mutate(All_surveys, Covid_ever = ifelse(rowSums(All_surveys[, 2:12], na.rm = TRUE) > 0, 1, 0))
sum(All_surveys$Covid_ever == 1, na.rm = TRUE) # 20 confirm infection during complete survey period


##### -------- Relate COVID antibody concentrations to COVID survey responses -------- #####

All_surveys$Covid_vis_date <- Data_Covid$Date
All_surveys$Contam_anti <- Data_Covid$infected


Weeknrs  <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_weeknrs.csv")
Weekdata <- merge(All_surveys, Weeknrs, by = "Subj_ID")
# Date of baseline survey was set to 21-4-2020 for everyone for convenience (or is March 11, 2020 used? )
# Determine date of first "1" to Contam in survey responses 
# Calculate if any 1 was prior to blood sample - if yes Contam_anti should be yes

