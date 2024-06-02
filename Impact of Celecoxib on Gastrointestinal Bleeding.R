###### STEP 1- import the data ##############
EP748_laboriginal  <- "C:/Users/Bobie/Desktop/SPH Fall 2023/EP 748/EP748_lab.csv"
data1 <- read.csv(EP748_laboriginal)
head(data1)
###### STEP 2- examine the data #############
View(data1)


install.packages("tableone")
library(tableone)

CreateTableOne(data=data1)  
View(data1)  
# more informative names 


names(data1)[2] <- 'Celecoxib'
names(data1)[3] <- 'Male'
names(data1)[4] <- 'Age'
names(data1)[5] <-'GI bleed history'
names(data1)[6] <-'Current GI protective agent use'
names(data1)[7] <-'Warfarin'
names(data1)[8] <-'Number of medications used'
names(data1)[9] <-'Peptic ulcer history'
names(data1)[10] <-'Osteoarthritis'
names(data1)[11] <-'Rheumatoid arthritis'
names(data1)[12] <-'Past GI protective agent use'
names(data1)[13] <-'Hospitalization history'
names(data1)[14] <-'COPD'
names(data1)[15] <-'Steroid use'
names(data1)[16] <-'GI bleeding'

## Get variables names to separate out continuous from categorical
dput(names(data1))

## Vector of variables to summarize
allVars <- c("Celecoxib", "Male", "Age", "GI bleed history", "Current GI protective agent use", 
             "Warfarin", "Number of medications used", "Peptic ulcer history", 
             "Osteoarthritis", "Rheumatoid arthritis", "Past GI protective agent use", 
             "Hospitalization history", "COPD", "Steroid use")


## Vector of categorical variables that need transformation
catVars <-  c("Celecoxib", "Male", "GI bleed history", "Current GI protective agent use", 
              "Warfarin", "Peptic ulcer history", 
              "Osteoarthritis", "Rheumatoid arthritis", "Past GI protective agent use", 
              "Hospitalization history", "COPD", "Steroid use")


## Create table- overall cohort 
Cohort1_table1 <- CreateTableOne(vars=allVars, data=data1, factorVars = catVars)
Treatments_table1 <- CreateTableOne(vars=allVars, strata='Celecoxib', data=data1, factorVars = catVars)

## Create table- by treatment status
table_Cohort1_table1 <- print(Cohort1_table1)
table_Treatments_table1 <- print(Treatments_table1)  

## save both the tables as csv files
write.csv(table_Cohort1_table1 , file = "Table 1 cohort.csv")
write.csv(table_Treatments_table1, file = "Table 1 by treatment.csv")  

write.csv(data1,file ="Testcoxinsas.csv")

###### STEP 3- Crude measures between treatment and outcome, not accounting for any confounding #############
install.packages("epiDisplay")
library(epiDisplay)

## incidence rates stratified by treatment groups

events <- rowsum(data1$`GI bleeding`, data1$Celecoxib)
pyears <- rowsum((data1$survt)/365, data1$Celecoxib)
ci.poisson(events, pyears, alpha=.05)

View(events)
## unadjusted hazard ratio using Cox proportional hazards model 


install.packages("survival")
library(survival)

unadj_model <- coxph(Surv(data1$survt, data1$`GI bleeding`) ~ data1$Celecoxib, method = "breslow") 
summary(unadj_model)

adj_model <- coxph(Surv(data1$survt, data1$`GI bleeding`) ~ data1$Celecoxib + data1$Male + data1$Age +
                     data1$`GI bleed history` + data1$`Current GI protective agent use` +
                     data1$Warfarin + data1$`Number of medications used`, method = "breslow")  
summary(adj_model)

###### Confounding adjustment approach 2- Propensity score matching #############

## PSM diagnostic 1- PS overlap
install.packages("ggplot2")
library(ggplot2)

# calculate PS in the sample before matching

psmodel <- glm(data1$Celecoxib ~ data1$Male + data1$Age +
  data1$`GI bleed history` + data1$`Current GI protective agent use` +
                 data1$Warfarin + data1$`Number of medications used`, family = "binomial", data=data1)
data1$prop.scores <- psmodel[["fitted.values"]]

# plot overlap before matching

cbPalette <- c("blue", "pink")
ggplot(data1, aes(x=prop.scores, fill = factor(Celecoxib), y=after_stat(scaled)))+
  geom_density(alpha = 0.6)+
  scale_x_continuous("Propensity score", limits=c(0,1)) +
  scale_fill_manual(name="Treatment", values=cbPalette, breaks = c("0", "1"),
                    labels = c("Celecoxib", "NSAID"))+
  scale_y_continuous("Scaled density") +
  theme(axis.title.x = element_text(size = 12, vjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.key = element_rect(colour = NA))+ ggtitle("before PS matching")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=15))

install.packages("MatchIt")
library(MatchIt)

m.out <- matchit( data1$Celecoxib ~ data1$Male + data1$Age +
                    data1$`GI bleed history` + data1$`Current GI protective agent use` +
                    data1$Warfarin + data1$`Number of medications used`, data=data1,
                  method="nearest",
                  caliper=0.025,
                  replace=FALSE)

matched <- match.data(m.out) 
View(matched)

# plot overlap after matching

ggplot(matched, aes(x=distance, fill = factor(Celecoxib), y=..scaled..))+
  geom_density(alpha = 0.6)+
  scale_x_continuous("Propensity score", limits=c(0,1)) +
  scale_fill_manual(name="Treatment", values=cbPalette, breaks = c("0", "1"),
                    labels = c("Celecoxib", "NSAID"))+
  scale_y_continuous("Scaled density") +
  theme(axis.title.x = element_text(size = 12, vjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.key = element_rect(colour = NA))+ ggtitle("After PS matching")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=15))


## PSM diagnostic 2- table of patient characteristics by treatment in the matched sample

library(tableone)

## Vector of variables to summarize
allVars <- c("Male", "Age", "GI bleed history", "Current GI protective agent use",
             "Warfarin", "Number of medications used", "Peptic ulcer history",
             "Osteoarthritis", "Rheumatoid arthritis", "Past GI protective agent use",
             "Hospitalization history", "COPD", "Steroid use")

## Vector of categorical variables that need transformation
catVars <-  c("Male", "GI bleed history", "Current GI protective agent use",
              "Warfarin", "Peptic ulcer history",
              "Osteoarthritis", "Rheumatoid arthritis", "Past GI protective agent use",
              "Hospitalization history", "COPD", "Steroid use")


matched_treatments_t1 <- CreateTableOne(vars = allVars, strata = 'Celecoxib', data = matched, factorVars = catVars)

table_matched_treatments_t1 <- print(matched_treatments_t1)

write.csv(table_matched_treatments_t1, file = "Table 1 PS matched.csv")

#analysis in the matched sample

PSM_model <- coxph(Surv(matched$survt, matched$`GI bleeding`)  ~ matched$Celecoxib, method = "breslow")  
summary(PSM_model)

###### Confounding adjustment approach 3- Propensity score weighting #############

data1$psw <- ifelse(data1$Celecoxib==1, 1/data1$prop.scores, 1/(1-data1$prop.scores))
View(data1)
#analysis in incorporating the IP weights

PSW_model <- coxph(Surv(data1$survt, data1$`GI bleeding`)  ~ data1$Celecoxib, weight= data1$psw, robust=TRUE, method = "breslow")  
summary(PSW_model)








