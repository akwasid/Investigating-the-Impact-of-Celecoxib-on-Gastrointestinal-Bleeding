 Define the file path and read the CSV file
EP748_laboriginal <- "C:/Users/Bobie/Desktop/SPH Fall 2023/EP 748/EP748_lab.csv"
data1 <- read.csv(EP748_laboriginal)
head(data1)
###### STEP 2 - Examine the data #############
View(data1)

# Install and load necessary packages
if (!require("tableone")) install.packages("tableone")
library(tableone)

# Create the initial table
CreateTableOne(data = data1)
View(data1)

# Rename columns for better readability
colnames(data1) <- c('ID', 'Celecoxib', 'Male', 'Age', 'GI_bleed_history', 
                     'Current_GI_protective_agent_use', 'Warfarin', 'Number_of_medications_used', 
                     'Peptic_ulcer_history', 'Osteoarthritis', 'Rheumatoid_arthritis', 
                     'Past_GI_protective_agent_use', 'Hospitalization_history', 'COPD', 
                     'Steroid_use', 'GI_bleeding', 'survt')

# Vector of variables to summarize
allVars <- c('Celecoxib', 'Male', 'Age', 'GI_bleed_history', 
             'Current_GI_protective_agent_use', 'Warfarin', 'Number_of_medications_used', 
             'Peptic_ulcer_history', 'Osteoarthritis', 'Rheumatoid_arthritis', 
             'Past_GI_protective_agent_use', 'Hospitalization_history', 'COPD', 
             'Steroid_use')

# Vector of categorical variables that need transformation
catVars <- c('Celecoxib', 'Male', 'GI_bleed_history', 'Current_GI_protective_agent_use', 
             'Warfarin', 'Peptic_ulcer_history', 'Osteoarthritis', 'Rheumatoid_arthritis', 
             'Past_GI_protective_agent_use', 'Hospitalization_history', 'COPD', 'Steroid_use')

# Create tables - overall cohort and by treatment status
Cohort1_table1 <- CreateTableOne(vars = allVars, data = data1, factorVars = catVars)
Treatments_table1 <- CreateTableOne(vars = allVars, strata = 'Celecoxib', data = data1, factorVars = catVars)

# Print and save tables as CSV files
write.csv(print(Cohort1_table1), file = "Table 1 cohort.csv")
write.csv(print(Treatments_table1), file = "Table 1 by treatment.csv")

write.csv(data1, file = "Testcoxinsas.csv")

###### STEP 3 - Crude measures between treatment and outcome, not accounting for any confounding #############

# Install and load epiDisplay package
if (!require("epiDisplay")) install.packages("epiDisplay")
library(epiDisplay)

# Incidence rates stratified by treatment groups
events <- rowsum(data1$GI_bleeding, data1$Celecoxib)
pyears <- rowsum(data1$survt / 365, data1$Celecoxib)
ci.poisson(events, pyears, alpha = .05)
View(events)

# Unadjusted hazard ratio using Cox proportional hazards model
if (!require("survival")) install.packages("survival")
library(survival)

unadj_model <- coxph(Surv(data1$survt, data1$GI_bleeding) ~ data1$Celecoxib, method = "breslow")
summary(unadj_model)

adj_model <- coxph(Surv(data1$survt, data1$GI_bleeding) ~ data1$Celecoxib + data1$Male + data1$Age +
                     data1$GI_bleed_history + data1$Current_GI_protective_agent_use +
                     data1$Warfarin + data1$Number_of_medications_used, method = "breslow")
summary(adj_model)

###### Confounding adjustment approach 2 - Propensity score matching #############

# Install and load ggplot2 and MatchIt packages
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

if (!require("MatchIt")) install.packages("MatchIt")
library(MatchIt)

# Calculate PS in the sample before matching
psmodel <- glm(Celecoxib ~ Male + Age + GI_bleed_history + Current_GI_protective_agent_use + 
                 Warfarin + Number_of_medications_used, family = "binomial", data = data1)
data1$prop_scores <- psmodel$fitted.values

# Plot overlap before matching
cbPalette <- c("blue", "pink")
ggplot(data1, aes(x = prop_scores, fill = factor(Celecoxib), y = after_stat(scaled))) +
  geom_density(alpha = 0.6) +
  scale_x_continuous("Propensity score", limits = c(0, 1)) +
  scale_fill_manual(name = "Treatment", values = cbPalette, breaks = c("0", "1"), labels = c("Celecoxib", "NSAID")) +
  scale_y_continuous("Scaled density") +
  theme(axis.title.x = element_text(size = 12, vjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
  theme(legend.key = element_rect(colour = NA)) +
  ggtitle("Before PS matching") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

# Perform matching
m.out <- matchit(Celecoxib ~ Male + Age + GI_bleed_history + Current_GI_protective_agent_use + 
                   Warfarin + Number_of_medications_used, data = data1, method = "nearest", caliper = 0.025, replace = FALSE)

matched <- match.data(m.out)
View(matched)

# Plot overlap after matching
ggplot(matched, aes(x = distance, fill = factor(Celecoxib), y = ..scaled..)) +
  geom_density(alpha = 0.6) +
  scale_x_continuous("Propensity score", limits = c(0, 1)) +
  scale_fill_manual(name = "Treatment", values = cbPalette, breaks = c("0", "1"), labels = c("Celecoxib", "NSAID")) +
  scale_y_continuous("Scaled density") +
  theme(axis.title.x = element_text(size = 12, vjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
  theme(legend.key = element_rect(colour = NA)) +
  ggtitle("After PS matching") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

# PSM diagnostic 2 - Table of patient characteristics by treatment in the matched sample
matched_treatments_t1 <- CreateTableOne(vars = allVars, strata = 'Celecoxib', data = matched, factorVars = catVars)
write.csv(print(matched_treatments_t1), file = "Table 1 PS matched.csv")

# Analysis in the matched sample
PSM_model <- coxph(Surv(matched$survt, matched$GI_bleeding) ~ matched$Celecoxib, method = "breslow")
summary(PSM_model)

###### Confounding adjustment approach 3 - Propensity score weighting #############

data1$psw <- ifelse(data1$Celecoxib == 1, 1 / data1$prop_scores, 1 / (1 - data1$prop_scores))
View(data1)

# Analysis incorporating the IP weights
PSW_model <- coxph(Surv(data1$survt, data1$GI_bleeding) ~ data1$Celecoxib, weights = data1$psw, robust = TRUE, method = "breslow")
summary(PSW_model)
