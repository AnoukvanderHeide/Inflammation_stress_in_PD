######################################################################
#### Perform linear regression with clinical variables, AvdH 2023 ####
######################################################################

# After QC, 594 samples are left and 66 proteins
# The limma package (linear models for microarray data), is best resource for differential expression analysis
# https://bioc.ism.ac.jp/packages/3.6/bioc/manuals/limma/man/limma.pdf

# To treat Subj_ID as random effect, use duplicateCorrelation function

# AveExpr: average log2-expression level for specific protein across all the arrays and channels
# Adj.p.value: p-value adjusted for multiple testing (Benjamini-Hochberg)
# logFC: log fold-change between groups

rm(list = ls())

#### -------- Prepare data -------- ####

# Load required packages
library(readr)
library(dplyr)
library(stats)
library(limma)
library(calibrate)
library(tidyverse)

# Load complete data frame with protein concentrations and clinical data
Data          <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink_Clinical_combined.csv")
IL8_index     <- which(names(Data) == "IL8")    # Get index of IL-8, which is the first protein in the dataset
Data          <- subset(Data, select = c(Subj_ID, SampleID, SR_class, Time, Visit, Age_cov, Sex, BMI, Smoking, DisDur, LEDD, STAI_t, BDI, QUIP, UPDRS2, UPDRS3_off, IL8_index:(ncol(Data) - 4)))

Data_Covid    <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_Covid_analysed.csv")   # file with COVID antibody results (including infected/vaccinated?)
Data          <- merge(Data, Data_Covid[, c("Subj_ID", "infected")], by = "Subj_ID", all.x = TRUE) # merge Olink data with the infected variable
Data$infected <- factor(Data$infected)          # Transform this variable into a factor
IL8_index     <- which(names(Data) == "IL8")    # Get index of IL-8, which is the first protein in the dataset
new_order     <- c(names(Data)[1:(IL8_index-1)], "infected", names(Data)[IL8_index:(ncol(Data)-1)]) 
Data          <- Data[, new_order]              # Add infected column before first protein
rm(Data_Covid, new_order, IL8_index)

# Handle missing covariate values
missing_age     <- Data$SampleID[is.na(Data$Age_cov)]  # No missings
missing_sex     <- Data$SampleID[is.na(Data$Sex)]      # No missings
missing_BMI     <- Data$SampleID[is.na(Data$BMI)]      # Missing for one patient for one visit, use the BMI of visit 1 for visit 2 as well
missing_smoking <- Data$SampleID[is.na(Data$Smoking)]  # No missings
missing_covid   <- Data$SampleID[is.na(Data$infected)] # Missings in 7 patients (15 samples)
Data$BMI[Data$SampleID == "sub-POMUA04706EA161294D4_vis2"] <- Data$BMI[Data$SampleID == "sub-POMUA04706EA161294D4_vis1"]
Data$infected[is.na(Data$infected)] <- "no"
rm(missing_age, missing_BMI, missing_smoking, missing_sex)

#### -------- Generate the design matrices -------- ####
# Define comparison of interest through a design matrix (covariates and grouping)

# Define factors and covariates
Age       <- Data$Age_cov
Sex       <- as.factor(Data$Sex)
BMI       <- Data$BMI
Infection <- Data$infected
#Smoking <- Data$Smoking # Not used now

Visits           <- as.factor(Data$Visit)
levels(Visits)   <- c("V1", "V2", "V3")
Covid            <- as.factor(Data$Time)
levels(Covid)    <- c("preCOVID", "inCOVID")
SR_class         <- as.factor(Data$SR_class)
levels(SR_class) <- c("low", "high")

# Design matrix for testing main effect of time (three visits separately)
Time                    <- factor(Visits, levels = c("V1", "V2", "V3"))
Design_visits           <- model.matrix(~ 0 + Time + Age + Sex + BMI + Infection)
colnames(Design_visits) <- c("V1", "V2", "V3", "Age", "Sex", "BMI", "Infection")
dim(Design_visits)      # Get the dimensions of design (594 7)

# Design matrix for testing main effect of time (pre vs. in-COVID)
Time                    <- factor(Covid, levels = c("preCOVID", "inCOVID"))
Design_cov              <- model.matrix(~ 0 + Time + Age + Sex + BMI + Infection)
colnames(Design_cov)    <- c("preCOVID", "inCOVID", "Age", "Sex", "BMI", "Infection")
dim(Design_cov)         # Get the dimensions of design (594 6)

# Design matrix for testing main effect of SR class (low vs. high)
SR                      <- factor(SR_class, levels = c("low", "high"))
Design_SR               <- model.matrix(~ 0 + SR + Age + Sex + BMI + Infection)
colnames(Design_SR)     <- c("low","high","Age", "Sex", "BMI", "Infection")
dim(Design_SR)          # Get the dimensions of design (594 6)

# Design matrix for interaction model
SR                      <- factor(SR_class, levels = c("low", "high"))
Time                    <- factor(Covid, levels = c("preCOVID", "inCOVID"))
Design_int              <- model.matrix(~ SR*Time + Age + Sex + BMI + Infection)
colnames(Design_int)    <- c("Intercept", "SR_class", "Time", "Age", "Sex", "BMI", "Infection", "Interaction")
dim(Design_int)         # Get the dimensions of design (594 8)


### -------- Generate the protein expression matrix -------- ###

# Transform NPX data for proteins to represent rows - no group names
IL8_index     <- which(names(Data) == "IL8")     # Get index of IL-8, which is the first protein in the dataset
Data_new      <- Data[, c(IL8_index:ncol(Data))] # Exclude columns with clinical information (only proteins remaining; starting at IL8)
dim(Data_new)                                    # Get the dimensions of new data set (594 66); first must match design dimensions
Data_new      <- as.matrix(t(Data_new))          # Transpose data
IDnr          <- as.factor(Data$Subj_ID)         # Number subjects
Subs          <- as.numeric(IDnr)                # Make vector of Subj_ID to use as random factor

# Fit the model for differentially expressed proteins for all three visits (independent of SR)
dupcor      <- duplicateCorrelation(Data_new, Design_visits, block=Subs) # To treat Subj_ID as random effect, use duplicateCorrelation function
fit         <- lmFit(Data_new, Design_visits, block=Subs, correlation=dupcor$consensus) # Fit the model
cont.matrix <- makeContrasts(V3-V2, V2-V1, V3-V1, levels = Design_visits) # Specify contrasts (what to compare with what)
fit_vis     <- contrasts.fit(fit,cont.matrix)                     # Estimate contrasts for each protein
fit_vis     <- eBayes(fit_vis)                                    # Calculate F-statistics
List_sig1   <- topTable(fit_vis, adjust = "BH", n = 15)           # List the top n most differentially expressed proteins
topTable(fit_vis,coef="V3 - V2") # MMP10 sign.
topTable(fit_vis,coef="V2 - V1") # Flt3L, uPA, CXCL1, TRAIL, ST1A1 and CXCL5 sign.
topTable(fit_vis,coef="V3 - V1") # MMP10, MMP1 sign.
summary(decideTests(fit_vis))

# Fit the model for differentially expressed proteins for pre vs in COVID (independent of SR)
dupcor       <- duplicateCorrelation(Data_new, Design_cov, block=Subs) # To treat Subj_ID as random effect, use duplicateCorrelation function
fit          <- lmFit(Data_new, Design_cov, block=Subs, correlation=dupcor$consensus) # Fit the model
cont.matrix  <- makeContrasts(inCOVID-preCOVID, levels = Design_cov)  # Specify contrasts (what to compare with what)
fit_cov      <- contrasts.fit(fit, cont.matrix)                  # Estimate contrasts for each protein
fit_cov      <- eBayes(fit_cov)                                  # Calculate F-statistics
List_sig2    <- topTable(fit_cov, adjust = "BH", n = 15)         # List the top n most differentially expressed proteins
summary(decideTests(fit_cov))                                    # MMP10 and MMP1 downregulated, ST1A1 upregulated

# Fit the model for differentially expressed proteins for low and high SR (independent of time)
dupcor       <- duplicateCorrelation(Data_new, Design_SR, block=Subs) # To treat Subj_ID as random effect, use duplicateCorrelation function
fit          <- lmFit(Data_new, Design_SR, block=Subs, correlation=dupcor$consensus) # Fit the model
cont.matrix  <- makeContrasts(low-high, levels = Design_SR)      # Specify contrasts (what to compare with what)
fit_SR       <- contrasts.fit(fit, cont.matrix)                  # Estimate contrasts for each protein
fit_SR       <- eBayes(fit_SR)                                   # Calculate F-statistics
List_sig3    <- topTable(fit_SR, adjust = "BH", n = 15)          # List the top n most differentially expressed proteins
summary(decideTests(fit_SR))                                     # No sign. proteins

# Fit the model for differentially expressed proteins for SR class * time interaction (final model)
dupcor       <- duplicateCorrelation(Data_new, Design_int, block=Subs)             # To treat Subj_ID as random effect, use duplicateCorrelation function
fit          <- lmFit(Data_new, Design_int, block=Subs, correlation=dupcor$consensus) # Fit the model
cont.matrix  <- makeContrasts(Time = Time,
                              SR_class = SR_class,
                              Interaction = Interaction,
                              levels = Design_int)
fit_int      <- contrasts.fit(fit, cont.matrix)                  # Estimate contrasts for each protein
fit_int      <- eBayes(fit_int)                                  # Calculate F-statistics
topTable(fit_int, coef="SR_class",    adjust = "BH")             # List top 10 most differentially expressed proteins for main effect SR class
topTable(fit_int, coef="Time",        adjust = "BH")             # List top 10 most differentially expressed proteins for main effect Time
topTable(fit_int, coef="Interaction", adjust = "BH")             # List top 10 most differentially expressed proteins for interaction
summary(decideTests(fit_int))


#### -------- Plot results --------- ####

library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Display number of differential expressed genes
myCol <- brewer.pal(3, "Accent")
vennDiagram(decideTests(fit_cov), include= c("up","down"), counts.col=c("darkorchid3","deeppink3"), lwd = 2,
            cex = 1, circle.col=myCol, mar=c(0.1, 0.1, 0.1, 0.1))
vennDiagram(decideTests(fit_int, adjust.method = "none"), include= c("up","down"), counts.col=c("darkorchid3","deeppink3"), lwd = 2,
            cex = 1, circle.col=myCol, mar=c(0.1, 0.1, 0.1, 0.1)) # same plot for interaction model without correction for multiple testing

# Make volcano plot with adjusted p-value for interaction model
List_sig <- topTable(fit_int, coef = "Time") 
List_sig$Protein = rownames(List_sig)
List_sig <- subset(List_sig, P.Value < 0.05)

Volc_int <- ggplot(List_sig, aes(x = logFC, y = -log10(adj.P.Val), label = Protein)) +
                    geom_point(aes(color = adj.P.Val < 0.05), size = 2.5) +
                    geom_text_repel(data = List_sig) +
                    theme_bw() +
                    xlab("logFC") +
                    ylab("-log10(p-value)") +
                    geom_vline(xintercept = 0, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    geom_hline(yintercept = 1.3, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("p < 0.05", "adj. p < 0.05")) +
                    theme(axis.text = element_text(size = 13, colour = "black"), 
                          axis.title = element_text(size = 15, colour = "black"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 13))

# Make volcano plot with adjusted p-value for separate effect covid (pre vs in and vis3 vs. vis1)
List1 <- topTable(fit_cov)
List1$Protein = rownames(List1)
List1 <- subset(List1, P.Value < 0.05) # Only show proteins with original p-value <0.05
List2 <- topTable(fit_vis,coef="V3 - V1")
List2$Protein = rownames(List2)
List2 <- subset(List2, P.Value < 0.05) # Only show proteins with original p-value <0.05

Volc_cov <- ggplot(List1, aes(x = logFC, y = -log10(adj.P.Val), label = Protein)) +
                    geom_point(aes(color = adj.P.Val < 0.05), size = 2.5) +
                    geom_text_repel(data = List1) +
                    theme_bw() +
                    xlab("logFC") +
                    ylab("-log10(p-value)") +
                    geom_vline(xintercept = 0, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    geom_hline(yintercept = 1.3, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("p < 0.05", "adj. p < 0.05")) +
                    theme(axis.text = element_text(size = 13, colour = "black"), 
                          axis.title = element_text(size = 15, colour = "black"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 13))
Volc_vis <- ggplot(List2, aes(x = logFC, y = -log10(adj.P.Val), label = Protein)) +
                    geom_point(aes(color = adj.P.Val < 0.05), size = 2.5) +
                    geom_text_repel(data = List2) +
                    theme_bw() +
                    xlab("logFC") +
                    ylab("-log10(p-value)") +
                    geom_vline(xintercept = 0, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    geom_hline(yintercept = 1.3, col = "gray50", linetype = "dotted", linewidth = 1.3) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("p < 0.05", "adj. p < 0.05")) +
                    theme(axis.text = element_text(size = 13, colour = "black"), 
                          axis.title = element_text(size = 15, colour = "black"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 13))

# X-Axis: log2 fold change of protein expression between pre and in-COVID
# Positive values = higher expression, negative values = lower expression
# Y-Axis (-log10(P.Value)): high values = more statistically significant differences (smaller p-values)

ggsave("Volcanoplot_cov.svg", plot = Volc_int, width=5, height=5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )
ggsave("Volcanoplot_cov_sep.svg", plot = Volc_cov, width=5, height=5.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )


### Plot data for sign. markers as boxplots ###

# Boxplots; not separated per SR class
Data$Time <- factor(Data$Time, levels = c("pre", "post"))
Data$SR_class <- factor(Data$SR_class)

p1 <- ggplot(data = Data, aes(x = Time, y = MMP1)) + 
              geom_boxplot(fill = "blue4", alpha = 0.4) +
              geom_point(position = position_jitter(width = 0.2), alpha = 0.4) +
              scale_x_discrete(labels=c("pre-COVID","in-COVID")) +
              geom_point(alpha=0.4) +
              labs(x = "PPP visit", y = "MMP-1 (NPX value)") +
              theme_minimal()
MMP1 <- p1 + theme(axis.text = element_text(size = 16, colour = "black"), 
                   axis.title = element_text(size = 16, colour = "black"))

p2 <- ggplot(data = Data, aes(x = Time, y = ST1A1)) + 
              geom_boxplot(fill = "blue4", alpha = 0.4) +
              geom_point(position = position_jitter(width = 0.2), alpha = 0.4) +
              scale_x_discrete(labels=c("pre-COVID","in-COVID")) +
              geom_point(alpha=0.4) +
              labs(x = "PPP visit", y = "ST1A1 (NPX value)") +
              theme_minimal()
ST1A1 <- p2 + theme(axis.text = element_text(size = 16, colour = "black"), 
                    axis.title = element_text(size = 16, colour = "black"))

p3 <- ggplot(data = Data, aes(x = Time, y = MMP10)) + 
              geom_boxplot(fill = "blue4", alpha = 0.4) +
              geom_point(position = position_jitter(width = 0.2), alpha = 0.4) +
              scale_x_discrete(labels=c("pre-COVID","in-COVID")) +
              geom_point(alpha=0.4) +
              labs(x = "PPP visit", y = "MMP-10 (NPX value)") +
              theme_minimal()
MMP10 <- p3 + theme(axis.text = element_text(size = 16, colour = "black"), 
                    axis.title = element_text(size = 16, colour = "black"))



# Boxplots; separated per SR class
p4 <- ggplot(data = Data, aes(x = SR_class, y = MMP1, fill=interaction(SR_class, Time), alpha=Time)) + 
              geom_boxplot() +
              scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
              scale_alpha_manual(values=c(1, 0.4)) +
              geom_point(position=position_jitterdodge(),alpha=0.4) +
              scale_x_discrete(labels=c("Low SR","High SR")) +
              labs(x = "", y = "MMP-1 (NPX value)", color = "Time") +
              theme_minimal()
MMP1_SR <- p4 + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))

p5 <- ggplot(data = Data, aes(x = SR_class, y = ST1A1, fill=interaction(SR_class, Time), alpha=Time)) + 
              geom_boxplot() +
              scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
              scale_alpha_manual(values=c(1, 0.4)) +
              geom_point(position=position_jitterdodge(),alpha=0.4) +
              scale_x_discrete(labels=c("Low SR","High SR")) +
              labs(x = "", y = "ST1A1 (NPX value)", color = "Time") +
              theme_minimal()
ST1A1_SR <- p5 + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))

ggsave("Boxplot_MMP1_cov.svg", plot = MMP1, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )
ggsave("Boxplot_ST1A1_cov.svg", plot = ST1A1, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )
ggsave("Boxplot_MMP10_cov.svg", plot = MMP10, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )
ggsave("Boxplot_MMP1_SR.svg", plot = MMP1_SR, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )
ggsave("Boxplot_ST1A1_SR.svg", plot = ST1A1_SR, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink plots" )





### Plots for Rick ###


p <- ggplot(data=Data, aes(x=SR_class, y=IL10, fill=interaction(SR_class, Time), alpha=Time))  +# resilience: 1=high, 2=low
  geom_boxplot() +
  scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
  scale_alpha_manual(values=c(1, 0.5, 0.2)) +
  geom_point(position=position_jitterdodge(),alpha=0.4) +          
  scale_x_discrete(labels=c("Low SR","High SR")) +
  xlab("") +
  ylab("IL6") +
  theme_minimal()
IL6 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))

p2 <- ggplot(data=Data, aes(x=SR_class, y=TNF, fill=interaction(SR_class, Time), alpha=Time))  +# resilience: 1=high, 2=low
  geom_boxplot() +
  scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
  scale_alpha_manual(values=c(1, 0.5, 0.2)) +
  geom_point(position=position_jitterdodge(),alpha=0.4) +          
  scale_x_discrete(labels=c("Low SR","High SR")) +
  xlab("") +
  ylab("TNF") +
  theme_minimal()
TNF <- p2 + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
