
rm(list = ls())

### -------- Prepare the data -------- ###

# Libraries
library(readr)
library(tidyr)
library(dplyr)
library(rstatix)
library(zoo)
library(nlme)
library(ggplot2)
library(ggpubr)

# Load all data files
Inflam   <- read_csv("M:/Documents/Projecten/Inflammatie/Analysis/Data_CRPNFL_wide_adj.csv")
PPP_data <- read_csv("M:/Documents/Projecten/Inflammatie/PPP_precovid_visit.csv")
Inflam   <- merge(Inflam, subset(PPP_data, select = c(Subj_ID, Age_cov, Sex, precov, DisDur, BMI, Smoking, LEDD, MoCA, BDI, STAI_s, QUIP, UPDRS3_off, HY_off)), by = "Subj_ID") # precov BDI an UPDRS3
Covid    <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/Covid_averages.csv") # to get SR and SR class

Inflam <- merge(Inflam, subset(Covid, select = c(Subj_ID, SR, SR_class)), by = "Subj_ID") # Add a column with continuous SR measure and a column with SR class (from latent class analysis; SR_spl2)
Inflam <- Inflam %>% rename(Age = Age_cov) # Rename age column
rm(Covid, PPP_data) # Remove data sets that are not needed anymore
  
# Identify missings in covariates
missing_BMI     <- Inflam$Subj_ID[is.na(Inflam$BMI)]     # no missings
missing_smoking <- Inflam$Subj_ID[is.na(Inflam$Sex)] # no missings
missing_age     <- Inflam$Subj_ID[is.na(Inflam$Age)]     # no missingS
missing_sex     <- Inflam$Subj_ID[is.na(Inflam$Sex)]     # no missings

# Make one data file in long format, with CRP, NFL and date in long format, while keeping the other variables from precovid visit
spec <- tribble(~ .name,    ~ .value, ~ timepoint,
                "CRP.Visit1", "CRP", "1",
                "CRP.Visit2", "CRP", "2",
                "CRP.Visit3", "CRP", "3",
                "NFL.Visit1", "NFL", "1",
                "NFL.Visit2", "NFL", "2",
                "NFL.Visit3", "NFL", "3",
                "Date.Visit1", "Date", "1",
                "Date.Visit2", "Date", "2",
                "Date.Visit3", "Date", "3",)

Inflam_long <- Inflam %>% pivot_longer_spec(spec, names_repair = "unique", values_drop_na = TRUE)
Inflam_long$Visit <- as.factor(Inflam_long$timepoint)
levels(Inflam_long$Visit) <- c("1", "2", "3")
Inflam_long$Sex <- as.factor(Inflam_long$Sex)
rm(spec)

# Show mean and SD for all 3 visits (whole group, no comparison resilience yet)
Inflam_long %>% group_by(Visit) %>% get_summary_stats(CRP, type = "mean_sd")
Inflam_long %>% group_by(Visit) %>% get_summary_stats(NFL, type = "mean_sd")

# Remove CRP values >20, then perform log transformation on marker concentrations
Inflam_long     <- subset(Inflam_long, CRP<20) 
Inflam_long$CRP <- log(Inflam_long$CRP)
Inflam_long <- Inflam_long %>% mutate(NFL_adj = 10^(log10(NFL) - 0.01103 * (Age - 18))) # Add age-adjusted NfL values
Inflam_long$NFL <- log(Inflam_long$NFL)

# Select group with 3 visits of which 2 pre-covid
Inflam_long$precov <- as.factor(Inflam_long$precov)
Inflam_long_pre2 <- Inflam_long %>% filter(precov == 2) # Select the group with 2 visits before covid
subj_count <- table(Inflam_long_pre2$Subj_ID) # Count occurrences of each Subj_ID
problematic <- names(subj_count[subj_count != 3]) # Check if there are any Subj_ID with fewer or more than three instances
print(problematic) # "sub-POMU1CF55A8FB405D10A" "sub-POMU782CA9DEB7DD68CC" "sub-POMUD69DA7AAD854ED06" "sub-POMUF8C5A88ADC866F24"

# Function to interpolate missing values of the 4 people with one data point missing
interpolate_within_group <- function(data) {
  data %>%
    group_by(Subj_ID) %>%
    mutate(CRP = ifelse(is.na(CRP), approx(Visit, CRP, method = "linear")$y, CRP),
           NFL = ifelse(is.na(NFL), approx(Visit, NFL, method = "linear")$y, NFL)) %>%
    ungroup()
}

Inflam_long_pre2_extended <- interpolate_within_group(Inflam_long_pre2) # Apply the interpolation function
new_rows <- Inflam_long_pre2 %>% # Extract data for Visit == 1 for problematic Subj_IDs
  filter(Visit == 1, Subj_ID %in% problematic)
new_rows$Visit <- as.factor(rep(2, nrow(new_rows))) # For all of these 4, visit 2 was missing
Inflam_long_pre2_extended <- bind_rows(Inflam_long_pre2_extended, new_rows) # Bind new rows with the extended dataframe
Inflam_long_pre2 <- arrange(Inflam_long_pre2_extended, Subj_ID, Visit) # Sort the dataframe
rm(new_rows, Inflam_long_pre2_extended)

# Prepare data frame for pre-post comparison
Inflam$CRP_1 <- ifelse(Inflam$precov==1, Inflam$CRP.Visit1, Inflam$CRP.Visit2) #For precov==1 use CRP of visit 1, voor precov==2 (all others) use CRP of visit 2
Inflam$CRP_2 <- ifelse(Inflam$precov==1, Inflam$CRP.Visit2, Inflam$CRP.Visit3) #For precov==1 use CRP of visit 2, voor precov==2 (all others) use CRP of visit 3
Inflam$NFL_1 <- ifelse(Inflam$precov==1, Inflam$NFL.Visit1, Inflam$NFL.Visit2) #For precov==1 use NFL of visit 1, voor precov==2 (all others) use NFL of visit 2
Inflam$NFL_2 <- ifelse(Inflam$precov==1, Inflam$NFL.Visit2, Inflam$NFL.Visit3) #For precov==1 use NFL of visit 2, voor precov==2 (all others) use NFL of visit 3
Inflam = Inflam[,!(names(Inflam) %in% c("CRP.Visit1","CRP.Visit2", "CRP.Visit3", # Remove remaining columns per visit, otherwise pivot_longer will not work properly
                                        "NFL.Visit1","NFL.Visit2", "NFL.Visit3",
                                        "Date.Visit1","Date.Visit2", "Date.Visit3"))]

# to wide format for selection pre/post visit
Inflam_wide <- Inflam %>% filter(!is.na(Inflam$CRP_1), !is.na(Inflam$CRP_2)) # Remove NA's
sum(is.na(Inflam_wide$CRP_1)) # Check if complete: should be 0
sum(is.na(Inflam_wide$CRP_2)) # Check if complete: should be 0

# Put back in long format
spec2 <- tribble(~ .name,    ~ .value, ~ timepoint,
                 "CRP_1", "CRP", "pre",
                 "CRP_2", "CRP", "post",
                 "NFL_1", "NFL", "pre",
                 "NFL_2", "NFL", "post")

Inflam_long_all       <- Inflam_wide %>% pivot_longer_spec(spec2, names_repair = "unique", values_drop_na = TRUE)
Inflam_long_all       <- Inflam_long_all %>% rename(Visit = timepoint)
Inflam_long_all$Visit <- as.factor(Inflam_long_all$Visit)
Inflam_long_all$CRP   <- as.double(Inflam_long_all$CRP)
Inflam_long_all$NFL   <- as.double(Inflam_long_all$NFL)
rm(spec2)

# Remove CRP values >20, then perform log transformation on marker concentrations
Inflam_long_all     <- subset(Inflam_long_all, CRP<20) 
Inflam_long_all$CRP <- log(Inflam_long_all$CRP)
#Inflam_long_all     <- Inflam_long_all %>% mutate(NFL_adj = 10^(log10(NFL) - 0.01103 * (Age - 18))) # Add age-adjusted NfL values
Inflam_long_all$NFL <- log(Inflam_long_all$NFL)

#### -------- Run linear mixed models for interaction stress-reactivity and COVID -------- ####

library(nlme)

### For group with 2 pre-covid visits ###

# Count number of patients included in analysis
length(unique(Inflam_long_pre2$Subj_ID)) # N=119

# Fit mixed-effects model with random intercept only
model_CRP_pre2 <- lme(CRP ~ SR_class * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long_pre2, na.action = na.omit)
model_NFL_pre2 <- lme(NFL ~ SR_class * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long_pre2, na.action = na.omit)

# Fit mixed-effects model with random slope for visit
model_CRP_pre2_slope <- lme(CRP ~ SR_class * Visit + Age + Sex + BMI, random= ~1 + Visit|Subj_ID, data=Inflam_long_pre2, na.action = na.omit)
model_NFL_pre2_slope <- lme(NFL ~ SR_class * Visit + Age + Sex + BMI, random= ~1 + Visit|Subj_ID, data=Inflam_long_pre2, na.action = na.omit)
#model_NFLadj_pre2_slope <- lme(NFL_adj ~ SR_class * Visit + Sex + BMI, random=~1 + Visit|Subj_ID, data=Inflam_long_pre2, na.action = na.omit)
AIC(model_CRP_pre2, model_CRP_pre2_slope) # Random slope model better
AIC(model_NFL_pre2, model_NFL_pre2_slope) # Random slope model better

anova(model_CRP_pre2_slope) # interaction p=.029
anova(model_NFL_pre2_slope) # nothing sign.
#anova(model_NFLadj_pre2_slope) # nothing sign.

# posthoc test for CRP
posthoc1 <- Inflam_long_pre2 %>% # post hoc test for SR class separately for each visit
            group_by(Visit) %>% 
            pairwise_t_test(CRP ~ SR_class, p.adjust.method="bonferroni")
posthoc2 <- Inflam_long_pre2 %>% # post hoc test for visit separately per SR class
            group_by(SR_class) %>%
            anova_test(dv = CRP, wid = Subj_ID, within = Visit) %>%
            get_anova_table() %>%
            adjust_pvalue(method = "bonferroni")
posthoc1 # Pairwise comparisons show that the effect SR class is only significant for visit 3
posthoc2 

# Function to calculate t-values
calculate_t_value <- function(data, visit, group1, group2, formula) {
    subset_data <- data %>% filter(Visit == visit, SR_class %in% c(group1, group2))  # Filter based on Visit and SR_class
    t_test_result <- t.test(formula = formula, data = subset_data)
    return(t_test_result$statistic)}

# Get t-values for reporting posthoc test
t_values <- posthoc1 %>% rowwise() %>%
            mutate(t_value = calculate_t_value(data = Inflam_long_pre2, 
                                               visit = as.character(Visit), 
                                               group1 = group1, group2 = group2, 
                                               formula = CRP ~ SR_class))
print(t_values) # t=-2.47 for visit 3; df=119-2 groups=117


### For whole group ###

# Count number of patients included in analysis
length(unique(Inflam_long_all$Subj_ID)) # N=214

# Fit mixed-effects model with random intercept only
model_CRP <- lme(CRP ~ SR_class * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long_all, na.action = na.omit)
model_NFL <- lme(NFL ~ SR_class * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long_all, na.action = na.omit)

# Fit mixed-effects model with random slope for visit
model_CRP_slope <- lme(CRP ~ SR_class * Visit + Age + Sex + BMI, random= ~1 + Visit|Subj_ID, data=Inflam_long_all, na.action = na.omit)
model_NFL_slope <- lme(NFL ~ SR_class * Visit + Age + Sex + BMI, random= ~1 + Visit|Subj_ID, data=Inflam_long_all, na.action = na.omit)
#model_NFLadj_slope <- lme(NFL_adj ~ SR_class * Visit + Sex + BMI, random= ~1 + Visit|Subj_ID, data=Inflam_long_all, na.action = na.omit)
AIC(model_CRP, model_CRP_slope) # Random slope model better
AIC(model_NFL, model_NFL_slope) # Random slope model better

anova(model_CRP_slope) # main effect SR_class sign. (p=.004) and interaction sign. (p=.033)
anova(model_NFL_slope) # nothing sign.
#anova(model_NFLadj_slope) # nothing sign.

# posthoc test for CRP
summary <-Inflam_long_all %>% group_by(SR_class, Visit) %>% get_summary_stats(CRP, type = "mean_sd")
data.frame(summary)
posthoc1 <- Inflam_long_all %>%      # post hoc test for SR class separately for each visit
            group_by(Visit) %>% 
            pairwise_t_test(CRP ~ SR_class, p.adjust.method="bonferroni")
posthoc2 <- Inflam_long_all %>%     # post hoc test for visit separately per SR class
            group_by(SR_class) %>%
            anova_test(dv = CRP, wid = Subj_ID, within = Visit) %>%
            get_anova_table() %>%
            adjust_pvalue(method = "bonferroni")
posthoc1 # Pairwise comparisons show that the effect SR class is only significant for post-COVID visit
posthoc2 # Not sign. for both groups

# Get t-values for reporting posthoc test
t_values <- posthoc1 %>% rowwise() %>% 
            mutate(t_value = calculate_t_value(data = Inflam_long_all, 
                                               visit = as.character(Visit), 
                                               group1 = group1, group2 = group2, 
                                               formula = CRP ~ SR_class))
print(t_values) #t(212)=-3.11, p=0.0012


#### -------- Run linear mixed models for relationship stress-symptoms -------- ####




# Mixed models for relationship between marker concentrations and clinical scores (including effect of time)
anova(lme(CRP ~ BDI * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit)) # effect BDI on CRP sign. (p=0.0065)
anova(lme(NFL ~ BDI * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ QUIP * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit)) # effect QUIP on CRP sign. (p=0.0069)
anova(lme(NFL ~ QUIP * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ STAI_s * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(NFL ~ STAI_s * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ UPDRS3_off * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(NFL ~ UPDRS3_off * Visit + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))

# Mixed models for relationship between marker concentrations and clinical scores (NOT including effect of time)
anova(lme(CRP ~ BDI + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit)) # effect BDI on CRP sign. (p=0.0064)
anova(lme(NFL ~ BDI + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ QUIP + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit)) # effect QUIP on CRP sign. (p=0.0068)
anova(lme(NFL ~ QUIP + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ STAI_s + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(NFL ~ STAI_s + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(CRP ~ UPDRS3_off + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))
anova(lme(NFL ~ UPDRS3_off + Age + Sex + BMI, random=~1|Subj_ID, data=Inflam_long, na.action = na.omit))

### Make and save box-plots ###

# CRP: Boxplot for all 3 visits divided by SR class
p <- ggplot(data=Inflam_long_pre2, aes(x=SR_class, y=CRP, fill=interaction(SR_class, Visit), alpha=Visit)) + # resilience: 1=high, 2=low
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.5, 0.2)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +          
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            xlab("") +
            ylab("ln(C-reactive protein [CRP])") +
            theme_minimal()
CRP_pre2 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(CRP_pre2)

# NfL: Boxplot for all 3 visits divided by SR class
p <- ggplot(data=Inflam_long_pre2, aes(x=SR_class, y=NFL, fill=interaction(SR_class, Visit), alpha=Visit)) + # resilience: 1=high, 2=low
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.5, 0.2)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +          
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            xlab("") +
            ylab("ln(Neurofilament light chain [NfL])") +
            theme_minimal()
NFL_pre2 <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(NFL_pre2)

# CRP: Boxplot for pre vs. post COVID-concentrations divided by SR class
p <- ggplot(data=Inflam_long_all, aes(x=SR_class, y=CRP, fill=interaction(SR_class, Visit), alpha=Visit)) + 
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.4, 1, 0.4)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            xlab("") +
            ylab("ln(C-reactive protein [CRP])") +
            theme_minimal()
CRP_all <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(CRP_all)

# NfL: Boxplot for pre vs. post COVID-concentrations divided by SR class
p <- ggplot(data=Inflam_long_all, aes(x=SR_class, y=NFL, fill=interaction(SR_class, Visit), alpha=Visit)) + 
            geom_boxplot() +
            scale_fill_manual(values=c("#009988", "#CC3311", "#009988", "#CC3311")) +
            scale_alpha_manual(values=c(1, 0.4, 1, 0.4)) +
            geom_point(position=position_jitterdodge(),alpha=0.4) +
            scale_x_discrete(labels=c("\n \n Low SR","\n \n High SR")) +
            xlab("") +
            ylab("ln(Neurofilament light chain [NfL])") +
            theme_minimal()
NFL_all <- p + theme(legend.position = "none", axis.text = element_text(size = 17, colour = "black"), axis.title.y = element_text(size = 18, margin = margin(r=10)))
print(NFL_all)

# Plots for association with clinical scores
Inflam_long2 = Inflam_long[!Inflam_long$Subj_ID =="sub-POMU900F78E54F00A78A",] # remove extreme outlier BDI

BDI_CRP_obs <- ggplot(Inflam_long2, aes(x = BDI, y = CRPorig)) +
                      geom_point(color = "#7F7F7F") +
                      geom_smooth(method = "lm", linewidth = 1, color = "black", alpha = 0.2, fill = "blue") +
                      #stat_cor(method = "pearson", label.x = 24) +
                      labs(x = "Depressive symptoms (BDI)", y = "ln(C-reactive protein concentration)") +
                      theme_minimal() + 
                      theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                            axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                            axis.text = element_text(size = 13))

QUIP_CRP_obs <- ggplot(Inflam_long2, aes(x = QUIP, y = CRP)) +
                      geom_point(color = "#7F7F7F") +
                      geom_smooth(method = "lm", linewidth = 1, color = "black", alpha = 0.2, fill = "blue") +
                      #stat_cor(method = "pearson", label.x = 24) +
                      labs(x = "Impulsitivy symptoms (QUIP)", y = "ln(C-reactive protein concentration)") +
                      theme_minimal() + 
                      theme(axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 5)),
                            axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 5)),
                            axis.text = element_text(size = 13))



# Save plots
ggsave("CRP_pre2.svg", plot = CRP_pre2, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("NFL_pre2.svg", plot = NFL_pre2, width=6, height=4.5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )

ggsave("CRP_all.svg", plot = CRP_all, width=4.5, height=6, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("NFL_all.svg", plot = NFL_all, width=4.5, height=6, path = "M:/Documents/Projecten/Inflammatie/Analysis" )

ggsave("BDI_CRP.svg", plot = BDI_CRP_obs, width=5, height=5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )
ggsave("QUIP_CRP.svg", plot = QUIP_CRP_obs, width=5, height=5, path = "M:/Documents/Projecten/Inflammatie/Analysis" )


# Get characteristics table for pre-COVID visit
library(tableone)
Variables    <- c("Age", "Sex", "BMI", "DisDur", "UPDRS3_off", "HY_off", "LEDD", "MoCA", "BDI", "STAI_s")
Categorical  <- c("Sex", "HY_off") #Define categorical variables

table1 <- CreateTableOne(vars = Variables, data = Inflam, factorVars = Categorical)
table1

# Get table for first in-COVID visit
PPP_allvisits <- read_csv("M:/Documents/Projecten/Covid-19 survey/Data and analysis/PEP data/Data files/POM_allvisits.csv")
PPP_incov     <- PPP_allvisits[(PPP_allvisits$Visit == PPP_allvisits$precov+1),]
PPP_incov     <- merge(PPP_incov, Inflam, by = "Subj_ID")

Variables     <- c("BMI", "UPDRS3_off", "HY_off", "LEDD", "MoCA", "BDI", "STAI_s")
Categorical   <- c("HY_off") #Define categorical variables

table2        <- CreateTableOne(vars = Variables, data = PPP_incov, factorVars = Categorical)
table2

# Get CRP and NfL concentrations (not log-transformed) for pre- and in-COVID visit
Variables    <- c("CRP_1", "CRP_2", "NFL_1", "NFL_2")
Inflam       <- Inflam %>% mutate(CRP_1 = if_else(CRP_1 > 20, NA_real_, CRP_1),
                                  CRP_2 = if_else(CRP_2 > 20, NA_real_, CRP_2))

table3       <- CreateTableOne(vars = Variables, data = Inflam)
table3
