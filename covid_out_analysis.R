### Example R Program to Run Analyses in Main Text
### written by Thomas A Murray <murra484@umn.edu> and Jared Huling
### Created on 8/16/22
library(tidyverse); library(parallel); library(mice)

#Set your working directory to the location of the dataset
#setwd("C:\\Users\\Thomas\\Box\\MET-Covid_Open_DSMB_Reports\\data\\data_to_share_publicly") 

#Read in publicly available .csv dataset
itt.dat = read.csv("covid_out_data.csv") %>% as_tibble() %>% mutate(across(where(is.character),factor))
mitt.dat = subset(itt.dat,in_MITT_analysis=="Yes") #modified intention-to-treat dataset

################################
### Primary Outcome Analysis ### 
################################
### Prepare mITT analysis datasets for imputation
pasimp.mitt.dat = mitt.dat %>% select(id,outcome_1_by_d14,outcome_2_by_d14,outcome_3_by_d14,outcome_4_by_d14,outcome_7_by_d14,
                                      sex,age,bmi,days_first_symptoms_to_study_initiation,screen_covid_vaccine,variant_period,
                                      metformin, ivermectin, fluvoxamine,screen_covid_treatment,starts_with("race_demo"),starts_with("pd_comorbidities"))

### Multiply Impute Missing Values (the parlmice code is written to run in parallel on 10 compute cores; this code could be run serially using mice instead, but it's much slower)
pasimp <- parlmice(pasimp.mitt.dat[,-1], defaultMethod = rep("pmm",4), maxit=40, cluster.seed=123, n.core=10, n.imp.core=5)

### Define outcome variable
#Primary outcome: Hypoxemia or Supplemental O2, ED Visit, Hospitalization, Death
y_var = paste("I(",paste(paste0(paste0("outcome_",c(1:4,7),"_by_d14")," == 'Yes'"),collapse=" | "),")",sep="")
#Secondary outcome: ED Visit, Hospitalization, Death
#y_var = paste("I(",paste(paste0(paste0("outcome_",c(3:4,7),"_by_d14")," == 'Yes'"),collapse=" | "),")",sep="")
#Secondary outcome: Hospitalization, Death
#y_var = paste("I(",paste(paste0(paste0("outcome_",c(4,7),"_by_d14")," == 'Yes'"),collapse=" | "),")",sep="")

### Metformin Analysis
# Fit and summarize Logistic Regression Model
met.fits = with(pasimp, glm(formula(paste0(y_var,"~ screen_covid_vaccine + metformin + ivermectin + fluvoxamine")), family = binomial()))
summary(pool(met.fits), conf.int = TRUE, exponentiate = TRUE)[3,] # summarize coefficients on the odds ratio scale

### Ivermectin Analysis
iver.fits = with(pasimp, glm(formula(paste0(y_var,"~ screen_covid_vaccine + metformin + ivermectin")), family = binomial(), 
                             subset = (mitt.dat$in_ivermectin_analysis == "Yes"))) ## subset only to concurrent control data for Iver analysis
summary(pool(iver.fits), conf.int = TRUE, exponentiate = TRUE)[4,] # summarize coefficients on the odds ratio scale

### Fluvoxamine Analysis
fluv.fits = with(pasimp, glm(formula(paste0(y_var,"~ screen_covid_vaccine + metformin + fluvoxamine")), family = binomial(), 
                             subset = (mitt.dat$in_fluvoxamine_analysis == "Yes"))) ## subset only to concurrent control data for Fluvox analysis
summary(pool(fluv.fits), conf.int = TRUE, exponentiate = TRUE)[4,] # summarize coefficients on the odds ratio scale

### These estimates correspond to those in Table 2 of the main paper (exact matching may not occur due to Monte Carlo error in the multiple imputation process, and choice variables and participants used for multiple imputation for each drug)







#########################
### Symptoms Analysis ###
#########################
library(geeM)

### Code to calculate symptom scores
### We create the total symptom score using the following point system: 
#   1) For the 10 symptoms on the daily log: "none" = 0 pts, "mild" = 1 pt, "moderate" = 2 pts, "severe" = 3 pts.
#   2) For loss of taste or smell: "same as usual" = 0 pts, "less than usual" = 1 pt, "no taste/smell" = 2 pts.
#   3) For vomiting: 0 times = 0 pts, 1-2 times = 1 pt, 3-4 times = 2 pts, and >4 times = 3 pts.
#   4) For diarrhea: 0 times = 0 pts, 1-2 times = 1 pt, 3-4 times = 2 pts, and >4 times = 3 pts. 
#   The point totals for vomiting and diarrhea were determined using the following Guidance for Industry document: 
#   "Assessing COVID-19-Related Symptoms in Outpatient Adult and Adolescent Subjects in Clinical Trials of Drugs and Biological Products for COVID-19 Prevention or Treatment Guidance for Industry."
#   We then computed the total symptom score for each participant by adding together the points corresponding to each of the above symptoms.

### transform to tall format
symptom_data_imp_tall <- mitt.dat %>%
  pivot_longer(
    cols = dplyr::contains("_symp_"), 
    names_to = c("day",".value"), 
    names_sep = "_symp_"
  ) %>%
  dplyr::mutate(day = as.numeric(gsub("^d", "", day))) %>%
  arrange(id, day)

### Determine points for each symptom on the daily symptom log and add to data set
symp_name <- c("Stuffy_or_runny_nose","Sore_throat","Shortness_of_breath_or_difficulty_breathing","Cough","Low_energy_or_tiredness","Muscle_or_body_ache","Headache","Chills_or_shivering","Feeling_hot_or_feverish","Nausea"); ten_symp <- c()
for(i in 1:10){
  foo <- 0*(symptom_data_imp_tall[,symp_name[i]]=="None") + 1*(symptom_data_imp_tall[,symp_name[i]]=="Mild") + 2*(symptom_data_imp_tall[,symp_name[i]]=="Moderate") + 3*(symptom_data_imp_tall[,symp_name[i]]=="Severe")
  ten_symp <- cbind(ten_symp, foo)
}
ten_symp <- setNames(data.frame(ten_symp), paste(symp_name,"pt",sep="_"))
symptom_data_imp_tall <- cbind(symptom_data_imp_tall, ten_symp)

### diarrhea, vomit are currently characters -> make numeric
symptom_data_imp_tall$diarrhea_number <- as.numeric(symptom_data_imp_tall$diarrhea_number)
symptom_data_imp_tall$vomit_number <- as.numeric(symptom_data_imp_tall$vomit_number)

### Determine points for loss of taste, loss of smell, diarrhea, vomit
symptom_data_imp_tall$vomit_pt <- 0*(symptom_data_imp_tall$vomit_number==0) + 1*(symptom_data_imp_tall$vomit_number>0) + 1*(symptom_data_imp_tall$vomit_number>2) + 1*(symptom_data_imp_tall$vomit_number>4)
symptom_data_imp_tall$diarrhea_pt <- 0*(symptom_data_imp_tall$diarrhea_number==0) + 1*(symptom_data_imp_tall$diarrhea_number>0) + 1*(symptom_data_imp_tall$diarrhea_number>2) + 1*(symptom_data_imp_tall$diarrhea_number>4)
symptom_data_imp_tall$smell_pt <- 0*(symptom_data_imp_tall$smell=="Same as usual") + 1*(symptom_data_imp_tall$smell=="Less than usual") + 2*(symptom_data_imp_tall$smell=="No taste/smell")
symptom_data_imp_tall$taste_pt <- 0*(symptom_data_imp_tall$taste_smell=="Same as usual") + 1*(symptom_data_imp_tall$taste_smell=="Less than usual") + 2*(symptom_data_imp_tall$taste_smell=="No taste/smell")

### Compute symptom scores
symptom_data_imp_tall$covid_specific_symptom_score <- rowSums(symptom_data_imp_tall[,grep("pt", colnames(symptom_data_imp_tall))][,-c(1,2,11,12,13,14,15)])
symptom_data_imp_tall$total_symptom_score <- rowSums(symptom_data_imp_tall[,grep("pt", colnames(symptom_data_imp_tall))][,-c(1,16)])
hist(symptom_data_imp_tall$total_symptom_score) # very right-skewed
hist(symptom_data_imp_tall$covid_specific_symptom_score) # very right-skewed
summary(symptom_data_imp_tall$total_symptom_score) 
summary(symptom_data_imp_tall$covid_specific_symptom_score) 



### Create dataset with day 1 scores as baseline score for adjustment
vars = c("id","total_symptom_score","covid_specific_symptom_score","day","screen_covid_vaccine","metformin","ivermectin","fluvoxamine","in_ivermectin_analysis","in_fluvoxamine_analysis")
blsymp = select(subset(symptom_data_imp_tall,day==1),c("id","total_symptom_score","covid_specific_symptom_score")) %>% mutate(bltot=total_symptom_score-mean(total_symptom_score,na.rm=TRUE),blcsp=covid_specific_symptom_score-mean(covid_specific_symptom_score,na.rm=TRUE),total_symptom_score=NULL,covid_specific_symptom_score=NULL)
fusymp = select(subset(symptom_data_imp_tall,day>1),all_of(vars))
symp.dat = left_join(fusymp,blsymp,by="id")

### Fit GEE Models
# Metformin Total Symptom Score (All Symptoms)
met.tot.fit = geem(total_symptom_score ~ factor(day) + bltot + screen_covid_vaccine + metformin + ivermectin + fluvoxamine, id = id, data = na.omit(symp.dat), family = gaussian(link = "identity"), corstr = "ar1")
summary(met.tot.fit) #Metformin Main Effect
met.tot.fit.byday = geem(total_symptom_score ~ factor(day)*metformin + bltot + screen_covid_vaccine + ivermectin + fluvoxamine, id = id, data = na.omit(symp.dat), family = gaussian(link = "identity"), corstr = "ar1")
summary(met.tot.fit.byday) #Metformin Effects by Day
# Metformin Covid-Specific Symptom Score 
met.covid.fit = geem(covid_specific_symptom_score ~ factor(day) + blcsp + screen_covid_vaccine + metformin + ivermectin + fluvoxamine, id = id, data = na.omit(symp.dat), family = gaussian(link = "identity"), corstr = "ar1")
summary(met.covid.fit) #Metformin Main Effect
met.covid.fit.byday = geem(covid_specific_symptom_score ~ factor(day)*metformin + blcsp + screen_covid_vaccine + ivermectin + fluvoxamine, id = id, data = na.omit(symp.dat), family = gaussian(link = "identity"), corstr = "ar1")
summary(met.tot.fit.byday) #Metformin Effects by Day

# Ivermectin Total Symptom Score (All Symptoms)
iver.tot.fit = geem(total_symptom_score ~ factor(day) + bltot + screen_covid_vaccine + metformin + ivermectin, id = id, data = na.omit(subset(symp.dat,in_ivermectin_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(iver.tot.fit) #Ivermectin Main Effect
iver.tot.fit.byday = geem(total_symptom_score ~ factor(day)*ivermectin + bltot + screen_covid_vaccine + metformin, id = id, data = na.omit(subset(symp.dat,in_ivermectin_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(iver.tot.fit.byday) #Ivermectin Effects by Day
# Ivermectin Covid-Specific Symptom Score 
iver.covid.fit = geem(covid_specific_symptom_score ~ factor(day) + blcsp + screen_covid_vaccine + metformin + ivermectin, id = id, data = na.omit(subset(symp.dat,in_ivermectin_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(iver.covid.fit) #Ivermectin Main Effect
iver.covid.fit.byday = geem(covid_specific_symptom_score ~ factor(day)*ivermectin + blcsp + screen_covid_vaccine + metformin, id = id, data = na.omit(subset(symp.dat,in_ivermectin_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(iver.covid.fit.byday) #Ivermectin Effects by Day

# Fluvoxamine Total Symptom Score (All Symptoms)
fluv.tot.fit = geem(total_symptom_score ~ factor(day) + bltot + screen_covid_vaccine + metformin + fluvoxamine, id = id, data = na.omit(subset(symp.dat,in_fluvoxamine_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(fluv.tot.fit) #Fluvoxamine Main Effect
fluv.tot.fit.byday = geem(total_symptom_score ~ factor(day)*fluvoxamine + bltot + screen_covid_vaccine + metformin, id = id, data = na.omit(subset(symp.dat,in_fluvoxamine_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(fluv.tot.fit.byday) #Fluvoxamine Effects by Day
# Fluvoxamine Covid-Specific Symptom Score 
fluv.covid.fit = geem(covid_specific_symptom_score ~ factor(day) + blcsp + screen_covid_vaccine + metformin + fluvoxamine, id = id, data = na.omit(subset(symp.dat,in_fluvoxamine_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(fluv.covid.fit) #Fluvoxamine Main Effect
fluv.covid.fit.byday = geem(covid_specific_symptom_score ~ factor(day)*fluvoxamine + blcsp + screen_covid_vaccine + metformin, id = id, data = na.omit(subset(symp.dat,in_fluvoxamine_analysis=="Yes")), family = gaussian(link = "identity"), corstr = "ar1")
summary(fluv.covid.fit.byday) #Fluvoxamine Effects by Day

### These GEE fits correspond to Figure 2 in the main paper
