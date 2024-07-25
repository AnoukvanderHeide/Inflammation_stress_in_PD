#### Exploratory analysis of Olink data in PPP, AvdH 2023

# This script removes samples with warnings and outliers, and markers with bad detection. Then the data is combined with
# relevant clinical variables and some exploratory statistical tests are performed.

# This script requires the Olink data in the exact format as in the NPX file. For this, use the adjusted xlsx file in 
# which the first and last rows with standard NPX info is copied directly from the original file. The rest of the data 
# is added from the PEP downloaded files, created in Structure_downloaded_PEP_files.R.

rm(list = ls())

library(readr)
library(stats)
library(dplyr)
library(ggplot2)
library(OlinkAnalyze)
library(stringr)
library(tidyr)
library(readxl)

# Load the adjusted file (which corresponds to the Routput file but than with the first and last rows from the original NPX file)
Data_olink <- read_NPX("M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink96_Adjusted.xlsx")

# PCA plot for 'QC warning' with each warning sample in red
Data_olink %>% 
    filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
    olink_pca_plot(df = ., color_g = "QC_Warning", byPanel = TRUE)  # generates PCA projection of all samples along two principal components

# Extract a list of all samples with QC-warnings to later exclude
warnings <- unique(Data_olink$SampleID[Data_olink$QC_Warning == "Warning"])

# PCA plot with outlier samples identified
QC_plot <- olink_pca_plot(df=Data_olink, color_g = "QC_Warning",
                          outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, quiet = TRUE)

# Generate box plots of NPX values for each sample, as initial QC step to identify potential warnings
Data_olink %>% 
    olink_dist_plot() +
    theme(axis.text.x = element_blank())

# Extract a list of outlier samples from the plot, to later exclude
outliers <- data.frame()
for (i in seq_along(QC_plot)) {plot_data <- QC_plot[[i]]$data
                               outliers <- rbind(outliers, 
                                                 plot_data[plot_data$Outlier == 1, c("SampleID", "Outlier", "Panel")])}
outliers <- unique(outliers)
print(outliers)

# Identify markers with >70% below the limit of detection
LOD_assays <- unique(Data_olink$Assay[Data_olink$MissingFreq > 0.70]) # Get list of markers for which >70% of the samples are <LOD
LOD_assays <- gsub("[- ]", "", LOD_assays)                            # Remove hyphens and spaces in marker names in NPX file
LOD_assays <- gsub("IL22RA1", "IL22.RA1", LOD_assays)                 # Change this name specifically, because for some reason a point is included in the name

# Load file with Olink data that was obtained in Structure_downloaded_PEP_files.R
Data_olink_orig <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink96_Routput.csv")
Data_olink_orig$SampleID <- Data_olink_orig$Subj_ID

# Load file with relevant clinical variables (in wide format)
Data_clinical <- read_csv("M:/Documents/Projecten/Inflammatie/POM_allvisits.csv")
Data_clinical$SampleID <- ifelse(Data_clinical$Visit == 1, paste(Data_clinical$Subj_ID, "_vis1", sep = ""),
                               ifelse(Data_clinical$Visit == 2, paste(Data_clinical$Subj_ID, "_vis2", sep = ""),
                                      ifelse(Data_clinical$Visit == 3, paste(Data_clinical$Subj_ID, "_vis3", sep = ""),
                                             Data_clinical$Subj_ID)))
Data_clinical$Time <- as.factor(ifelse(Data_clinical$precov == 1 & Data_clinical$Visit == 1, "pre",
                                       ifelse(Data_clinical$precov == 2 & (Data_clinical$Visit == 1 | Data_clinical$Visit == 2), "pre",
                                              ifelse(Data_clinical$precov == 3 & Data_clinical$Visit == 3, "pre",
                                                     "post"))))
Data_clinical$SR_class <- as.factor(Data_clinical$SR_class)

# Combine clean Olink markers data file (in wide format) with relevant clinical data
Data_combined_wide <- merge(Data_clinical[, -which(names(Data_clinical) == "Visit_date")], Data_olink_orig[, -which(names(Data_olink_orig) == "Subj_ID")], by="SampleID") # 616 samples

# Remove QC warning and outlier samples, remove markers 
Data_combined_wide <- Data_combined_wide %>% filter(!SampleID %in% warnings) # Remove warnings: 598 samples left (18 removed)
Data_combined_wide <- Data_combined_wide %>% filter(!SampleID %in% outliers$SampleID) # Remove outliers: 594 samples left (4 additional samples removed)
Data_combined_wide <- Data_combined_wide[, !colnames(Data_combined_wide) %in% LOD_assays] # Remove markers with >70% below LOD: 26 must be removed

# To add exact date of sample, use CRPNFL file instead of clinical data file (which is based on week number and is less specific)
Dates <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_CRPNFL_long.csv")[ ,c("Subj_ID", "Visit", "Date")]
Dates$SampleID <- ifelse(Dates$Visit == 1, paste(Dates$Subj_ID, "_vis1", sep = ""),
                                 ifelse(Dates$Visit == 2, paste(Dates$Subj_ID, "_vis2", sep = ""),
                                        ifelse(Dates$Visit == 3, paste(Dates$Subj_ID, "_vis3", sep = ""),
                                               Dates$Subj_ID)))
Data_combined_wide <- merge(Data_combined_wide, Dates[, c("SampleID", "Date")], by="SampleID")
Data_combined_wide <- Data_combined_wide %>% relocate(Date, .before = SR_class)

#write.csv(Data_combined_wide, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink_Clinical_combined.csv", row.names = FALSE)

# Create long format as well for exploratory tests
Data_olink_filtered <- Data_olink %>% filter(!SampleID %in% warnings) # Remove warnings: 598 samples left
Data_olink_filtered <- Data_olink %>% filter(!SampleID %in% outliers$SampleID)
Data_olink_filtered <- Data_olink[Data_olink$MissingFreq <= 0.70, ]
Data_combined_long  <- merge(Data_clinical[, c("Subj_ID", "SampleID", "Time", "Visit", "precov", "SR_class", "Age_cov", "Sex", "BMI", "Smoking", "BDI", 
                                               "STAI_t", "QUIP", "DisDur", "LEDD", "UPDRS2", "UPDRS3_off")], Data_olink_filtered, by="SampleID")

write.csv(Data_combined_long, file = "M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink_Clinical_combined_long.csv", row.names = FALSE)

### -------- Exploratory statistical tests -------- ###

# Compare pre-covid and in-covid visit samples: no markers with sign difference
ttest_time  <- olink_ttest(df = Data_combined_long, variable = "Time", alternative = "two.sided")
top_10_time <- ttest_time %>% slice_head(n = 10) %>% pull(OlinkID) # select names of the top 10 most significant proteins
olink_volcano_plot(p.val_tbl = ttest_time,                         # volcano plot with annotated top 10 most significant proteins
                   x_lab = 'In-covid vs. pre-covid',
                   olinkid_list = top_10_time)

# Compare high and low SR samples: no markers with sign difference
ttest_SR    <- olink_ttest(df = Data_combined_long, variable = "SR_class", alternative = "two.sided")
top_10_SR   <- ttest_SR %>% slice_head(n = 10) %>% pull(OlinkID) # select names of the top 10 most significant proteins
olink_volcano_plot(p.val_tbl = ttest_SR,                         # volcano plot with annotated top 10 most significant proteins
                   x_lab = 'High vs. low SR class',
                   olinkid_list = top_10_SR)


### -------- PCA plots separate per stressor-reactivity class -------- ### 

PCA_highSR <- Data_combined_wide %>% filter(!grepl("1", SR_class, ignore.case = TRUE)) # all people with SR 2 
PCA_lowSR <- Data_combined_wide %>% filter(!grepl("2", SR_class, ignore.case = TRUE))  # all people with SR 1

  # High SR group: Plot all 3 visits together in a scatter plot with elipses
  library(cowplot)
  
  PCA_highSR$Time <- as.numeric(PCA_highSR$Time) # kan niet want pre en post niet 1 en 2
  markers <- PCA_highSR[35:100]  #34-125 are markers, column 33 is Time and column 1 is SampleID
  pc <- prcomp(markers)
  df <- cbind(pc$x, PCA_highSR[,1], PCA_highSR[,5]) %>% as.data.frame()
  df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(PCA_highSR)))
  df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(PCA_highSR)))
  df$V68 <- as.factor(df$V68)  # v68 is the name of the column with visit number (check in df)
  
  p1 <- ggplot(df, aes(PC1, PC2, colour = V68)) +
               geom_point(size = 3, aes(shape = V68)) +
               stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
                            data = df[df$V68 == "1" | df$V68 == "2" | df$V68 == "3",], size = 1)
  
  # Add density curves to the x and y axes
  xdens <- axis_canvas(p1, axis = "x") + 
                       geom_density(data = df, aes(x = PC1, fill = V68, colour = V68), alpha = 0.3)
  ydens <- axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
                       geom_density(data = df, aes(x = PC2, fill = V68, colour = V68), alpha = 0.3) +
                       coord_flip()
  p1 %>% insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
         insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
         ggdraw()
  
  # Exploratory regression with factors time (pre vs. in-covid) and stress-reactivity class (high vs. low)
  model <- olink_lmer(df = Data_combined_long, 
                                    variable = c('Time', 'SR_class'),
                                    random = 'Time')
  
  # Extracting the significant proteins: nothing significant
  model_sig <- model %>%
    filter(Threshold == 'Significant', term == 'SR_class') %>%
    pull(OlinkID)
  
  # Low SR group: Plot all 3 visits together in a scatter plot with elipses
  PCA_lowSR$Time <- as.numeric(PCA_lowSR$Time) # kan niet want pre en post niet 1 en 2
  markers <- PCA_lowSR[35:100]  #34-125 are markers, column 33 is Time and column 1 is SampleID
  pc <- prcomp(markers)
  df <- cbind(pc$x, PCA_lowSR[,1], PCA_lowSR[,5]) %>% as.data.frame()
  df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(PCA_lowSR)))
  df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(PCA_lowSR)))
  df$V68 <- as.factor(df$V68)  # v68 is the name of the column with visit number (check in df)
  
  p2 <- ggplot(df, aes(PC1, PC2, colour = V68)) +
    geom_point(size = 3, aes(shape = V68)) +
    stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
                 data = df[df$V68 == "1" | df$V68 == "2" | df$V68 == "3",], size = 1)
  
  # Add density curves to the x and y axes
  xdens <- axis_canvas(p2, axis = "x") + 
                       geom_density(data = df, aes(x = PC1, fill = V68, colour = V68), alpha = 0.3)
  ydens <- axis_canvas(p2, axis = "y", coord_flip = TRUE) + 
                       geom_density(data = df, aes(x = PC2, fill = V68, colour = V68), alpha = 0.3) +
                       coord_flip()
  p2 %>% insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
         insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
         ggdraw()
