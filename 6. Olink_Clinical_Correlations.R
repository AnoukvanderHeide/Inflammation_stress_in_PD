########################################################################################
####    Calculate correlation between proteins and clinical variables, AvdH 2023    ####
########################################################################################

rm(list = ls())

# Load the required packages
library(readr)
library(stats)
library(dplyr)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(stringr)
library(corrr)
library(RcmdrMisc) # to use the function rcorr.adjust

# Load complete data frame with protein concentrations and clinical data
Data <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_Olink_Clinical_combined.csv")
IL8_index <- which(names(Data) == "IL8") 
cols_to_select <- c("Subj_ID", "SampleID", "SR_class", "Time", "Age_cov", "Sex", "BMI", "Smoking", "DisDur", "LEDD", "MoCA", "STAI_t", "BDI", "QUIP", "UPDRS2", "UPDRS3_off")
cols_to_select <- c(cols_to_select, names(Data)[IL8_index:(ncol(Data) - 4)]) # IL-8 is the first Olink value, select several clinical variables and then from IL-8 onwards (except last 4; QC warnings and plate nr)
Data <- Data[, cols_to_select]

#### -------- Prepare data -------- ####

# Handle missing data
missing_BMI <- Data$SampleID[is.na(Data$BMI)] # Missing for one patient for one visit, use the BMI of visit 1 for visit 2 as well
Data$BMI[Data$SampleID == "sub-POMUA04706EA161294D4_vis2"] <- Data$BMI[Data$SampleID == "sub-POMUA04706EA161294D4_vis1"]
missing_smoking <- Data$SampleID[is.na(Data$Smoking)] # no missings
missing_age <- Data$SampleID[is.na(Data$Age_cov)]     # no missingS
missing_sex <- Data$SampleID[is.na(Data$Sex)]         # no missings

# Mean imputation of clinical covariates
NAs <- sum(apply(Data, 1, function(row) any(is.na(row))))  # 50 rows have NAs
Data$STAI_t[is.na(Data$STAI_t)] <- mean(Data$STAI_t,na.rm=T)
Data$BDI[is.na(Data$BDI)] <- mean(Data$BDI,na.rm=T)
Data$QUIP[is.na(Data$QUIP)] <- mean(Data$QUIP,na.rm=T)
Data$LEDD[is.na(Data$LEDD)] <- mean(Data$LEDD,na.rm=T)
Data$UPDRS2[is.na(Data$UPDRS2)] <- mean(Data$UPDRS2,na.rm=T)
Data$UPDRS3_off[is.na(Data$UPDRS3_off)] <- mean(Data$UPDRS3_off,na.rm=T)
NAs <- sum(apply(Data, 1, function(row) any(is.na(row)))) # Now, 0 rows have NAs

# Make all variables numeric (needed for rcorr), so change factors time and sex to numbers
#Data$Sex <- as.numeric(ifelse(Data$Sex == "Male", 1, ifelse(Data$Sex == "Female", 2, Data$Sex)))
#Data$Time <- as.numeric(ifelse(Data$Time == "pre", 1, ifelse(Data$Time == "post", 2, Data$Time)))

Data_plot <- Data[,c(-1),] # exclude column with subject ID
#Data_plot <- Data[,c(-(3:5)),] # exclude column with BMI, age and sex (not interesting for final plots)

#### -------- Create subgroups based on time and SR -------- ####

# Divide the data frame based on time (pre vs. in-covid)
data_precov <- Data_plot %>% filter(!str_detect(Time, '2')) # pre-COVID visit: remove all where time=1
data_incov  <- Data_plot %>% filter(!str_detect(Time, '1'))
data_precov <- data_precov[,c(-1, -2),]     # Excluding column with time and SR_class
data_incov  <- data_incov[,c(-1, -2),]

# Divide the data frame based on SR_class (high vs. low)
data_lowSR <- Data_plot %>% filter(!str_detect(SR_class, '2')) # low SR: remove all where SR class=2
data_highSR  <- Data_plot %>% filter(!str_detect(SR_class, '1'))
data_lowSR <- data_lowSR[,c(-1, -2),]    # Excluding column with time and SR_class
data_highSR  <- data_highSR[,c(-1, -2, -3, -4, -5),]
#data_all    <- Data_plot[,c(-1, -2),]
data_all    <- Data_plot[,c(-(1:8)),]

### -------- Heat map correlations between clinical variables and proteins -------- ###

cor_precov <- rcorr(as.matrix(data_precov))
cor_incov  <- rcorr(as.matrix(data_incov))
cor_highSR  <- rcorr(as.matrix(data_highSR))
cor_lowSR <- rcorr(as.matrix(data_lowSR))
cor_all <- rcorr(as.matrix(data_all))

# Plot correlation heat maps per group (divided in time and SR class)
map1 <- corrplot(cor_precov$r[1:10, 11:ncol(data_precov)], # first argument is y axis (clinical measures) and second is x axis (proteins) in the new file
                 p.mat=cor_precov$P[1:10, 11:ncol(data_precov)], 
                 tl.col = "black", #color of the labels
                 method="color", 
                 col=colorRampPalette(c("turquoise4","white","lightcoral"))(100),
                 sig.level = c(0.001, 0.01, 0.05),
                 pch.cex = 0.9, #Size of the stars in the significant quadrants
                 insig = 'label_sig')

map2 <- corrplot(cor_incov$r[1:10, 11:ncol(data_incov)], # first argument is y axis (clinical measures) and second is x axis (proteins) in the new file
                 p.mat=cor_incov$P[1:10, 11:ncol(data_incov)], 
                 tl.col = "black", #color of the labels
                 method="color", 
                 col=colorRampPalette(c("turquoise4","white","lightcoral"))(100),
                 sig.level = c(0.001, 0.01, 0.05),
                 pch.cex = 0.9, #Size of the stars in the significant quadrants
                 insig = 'label_sig')

map3 <- corrplot(cor_highSR$r[1:10, 11:ncol(data_highSR)], # first argument is y axis (clinical measures) and second is x axis (proteins) in the new file
                 p.mat=cor_highSR$P[1:10, 11:ncol(data_highSR)], 
                 tl.col = "black", #color of the labels
                 method="color", 
                 col=colorRampPalette(c("turquoise4","white","lightcoral"))(100),
                 sig.level = c(0.001, 0.01, 0.05),
                 pch.cex = 0.9, #Size of the stars in the significant quadrants
                 insig = 'label_sig')

map4 <- corrplot(cor_lowSR$r[1:7, 8:ncol(data_lowSR)], # first argument is y axis (clinical measures) and second is x axis (proteins) in the new file
                 p.mat=cor_lowSR$P[1:7, 8:ncol(data_lowSR)], 
                 tl.col = "black", #color of the labels
                 method="color", 
                 col=colorRampPalette(c("turquoise4","white","lightcoral"))(100),
                 sig.level = c(0.001, 0.01, 0.05),
                 pch.cex = 0.9, #Size of the stars in the significant quadrants
                 insig = 'label_sig')

# Plot correlation heat maps for whole sample
map5 <- corrplot(cor_all$r[1:10, 11:ncol(data_all)], # first argument is y axis (clinical measures) and second is x axis (proteins) in the new file
                 p.mat=cor_all$P[1:10, 11:ncol(data_all)], 
                 tl.col = "black", #color of the labels
                 method="color", 
                 col=colorRampPalette(c("turquoise4","white","lightcoral"))(100),
                 sig.level = c(0.001, 0.01, 0.05),
                 pch.cex = 0.9, #Size of the stars in the significant quadrants
                 insig = 'label_sig')

# p<.01 BDI: IL6, IL18R1, CXCL5, CCL3, CXCL6, TSFRSF9, CSF1
# p<.01 STAI: MMP10, TNF, TSFRSF9

#### -------- Explore clinical correlations with partial correlations and mixed-effects models -------- ####
library(nlme)
library(lme4)
library(ggpubr)
library(ppcor)
library(ggeffects)

Data2 <- Data[Data$SampleID != "sub-POMU900F78E54F00A78A_vis2", ] # Extreme outlier in BDI, exclude from models

# Calculate partial correlations (taking covariates age sex and BMI into account)
cor_BDI_TNF <- pcor.test(Data2$BDI, Data2$TNF, Data2[, c("BMI", "Age_cov", "Sex")], method = "pearson")
print(cor_BDI_TNF$estimate) # R=0.13; p=.001
print(cor_BDI_TNF$p.value)
cor_BDI_IL6 <- pcor.test(Data2$BDI, Data2$IL6, Data2[, c("BMI", "Age_cov", "Sex")], method = "pearson")
print(cor_BDI_IL6$estimate) # R=0.09; p=.02
print(cor_BDI_IL6$p.value)

cor_STAI_TNF <- pcor.test(Data$STAI_t, Data$TNF, Data[, c("BMI", "Age_cov", "Sex")], method = "pearson")
print(cor_STAI_TNF$estimate) # R=0.18; p<.0001
print(cor_STAI_TNF$p.value)
cor_STAI_IL6 <- pcor.test(Data$STAI_t, Data$IL6, Data[, c("BMI", "Age_cov", "Sex")], method = "pearson")
print(cor_STAI_IL6$estimate) # R=0.08; p=.05
print(cor_STAI_IL6$p.value)

# Plot separate correlations between markers and psychological scores (not corrected and no covariates, just for visualization)
p <- ggscatter(Data, x = "STAI_t", y = "TNF", 
                add = "reg.line", conf.int = TRUE, 
                add.params = list(color = "#CC3311"),
                cor.coef = TRUE, cor.method = "pearson", # correlation 0.18, p<.0001
                cor.coeff.args = list(label.x = 60),
                xlab = "Anxiety symptoms (STAI)", ylab = "TNF-alpha (NPX value)")
STAI_TNF <- p + theme(axis.text = element_text(size = 17, colour = "black"), 
                 axis.title = element_text(size = 17, colour = "black"))

p <- ggscatter(Data, x = "STAI_t", y = "IL6", 
                add = "reg.line", conf.int = TRUE, 
                add.params = list(color = "#CC3311"),
                cor.coef = TRUE, cor.method = "pearson", # correlation 0.08, p=.05
                cor.coeff.args = list(label.x = 60),
                xlab = "Anxiety symptoms (STAI)", ylab = "IL-6 (NPX value)")
STAI_IL6 <- p + theme(axis.text = element_text(size = 17, colour = "black"), 
                      axis.title = element_text(size = 17, colour = "black"))

p <- ggscatter(Data2, x = "BDI", y = "TNF", 
                add = "reg.line", conf.int = TRUE, 
                add.params = list(color = "#CC3311"),
                cor.coef = TRUE, cor.method = "pearson", # correlation 0.13, p=.001
                cor.coeff.args = list(label.x = 24),
                xlab = "Depressive symptoms (BDI)", ylab = "TNF-alpha (NPX value)")
BDI_TNF <- p + theme(axis.text = element_text(size = 17, colour = "black"), 
                      axis.title = element_text(size = 17, colour = "black"))

ggscatter(Data2, x = "BDI", y = "IL6", 
                add = "reg.line", conf.int = TRUE, 
                add.params = list(color = "black"),  # Change line color to black
                conf.int.params = list(fill = "#CC3311"),  # Change CI color to #CCCCFF
                cor.coef = TRUE, cor.method = "pearson",
                cor.coeff.args = list(label.x = 24),
                xlab = "Depressive symptoms (BDI)", ylab = "Observed IL-6 concentration (NPX)") +
                theme_minimal() + 
                theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                      axis.text = element_text(size = 13))

ggscatter(Data, x = "QUIP", y = "TNF", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", # correlation 0.10, p=.016
          xlab = "Impulsivity (QUIP)", ylab = "NPX TNF-alpha")

ggsave("STAI_TNF.svg", plot = STAI_TNF, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("STAI_IL6.svg", plot = STAI_IL6, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("BDI_TNF.svg", plot = BDI_TNF, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("BDI_IL6.svg", plot = BDI_IL6, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )


# Run linear mixed models for relationship a-priori markers with BDI
TNF  <- anova(lme(TNF ~ BDI + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data2, na.action = na.omit))      # p=0.1314
IL6  <- anova(lme(IL6 ~ BDI + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data2, na.action = na.omit))      # p=0.0116
MCP1 <- anova(lme(MCP1 ~ BDI + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data2, na.action = na.omit))     # p=0.1866
IL10 <- anova(lme(IL10 ~ BDI + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data2, na.action = na.omit))     # p=0.3458

# Calculate adjusted p-values
pvalues<-c(0.0116, 0.1314, 0.1866, 0.3458) # p-values of effect of BDI on markers, ordered from low to high
fdrs<-p.adjust(pvalues, method="BH")
print(fdrs) # 0.0464 0.2488 0.2488 0.3458

# Linear mixed models for relationship a-priori markers with STAI
TNF_2 <- anova(lme(TNF ~ STAI_t + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data, na.action = na.omit))    # p=0.0507
IL6_2 <- anova(lme(IL6 ~ STAI_t + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data, na.action = na.omit))    # p=0.1630
MCP1_2 <- anova(lme(MCP1 ~ STAI_t + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data, na.action = na.omit))  # p=0.1308
IL10_2 <- anova(lme(IL10 ~ STAI_t + Age_cov + Sex + BMI, random=~1|Subj_ID, data=Data, na.action = na.omit))  # p=0.2473

# Calculate adjusted p-values
pvalues<-c(0.0507, 0.1308,  0.1630, 0.2473)
fdrs<-p.adjust(pvalues, method="BH", n = length(pvalues))
print(fdrs) # 0.203 0.217 0.217 0.247

# Plot predicted relationship IL-6 and BDI (only significant result mixed model)
Data2$Sex <- as.factor(Data2$Sex)
IL6_new <- lmer(IL6 ~ BDI + Age_cov + Sex + BMI + (1|Subj_ID), data=Data2, na.action = na.omit)  # Make same model IL6 with lme4, ggpredict does not work with nlme
IL6_pred <- ggpredict(IL6_new, terms = c("BDI"), condition = c(Sex = "Male"))                       # Use ggpredict to generate predictions
# ggpredict holds non-focal terms constant at mean value or at reference level (for factors). Here, we used the most common levels for sex (male), mean age and mean BMI

# Plot predicted IL6
BDI_IL6_pred <- ggplot(IL6_pred, aes(x = x, y = predicted)) +
                       geom_line(linewidth = 1) +
                       geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "blue") +
                       labs(x = "Depressive symptoms (BDI)", y = "Predicted IL-6 concentration (NPX)") +
                       theme_minimal() + 
                       theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                             axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                             axis.text = element_text(size = 13))

# Plot observed IL6
BDI_IL6_obs <- ggplot(Data2, aes(x = BDI, y = IL6)) +
                      geom_point(color = "#7F7F7F") +
                      geom_smooth(method = "lm", linewidth = 1, color = "black", alpha = 0.2, fill = "blue") +
                      stat_cor(method = "pearson", label.x = 24) +
                      labs(x = "Depressive symptoms (BDI)", y = "Observed IL-6 concentration (NPX)") +
                      theme_minimal() + 
                      theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                            axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                            axis.text = element_text(size = 13))

ggsave("BDI_IL6_pred.svg", plot = BDI_IL6_pred, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink" )
ggsave("BDI_IL6_obs.svg", plot = BDI_IL6_obs, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis/Olink" )
