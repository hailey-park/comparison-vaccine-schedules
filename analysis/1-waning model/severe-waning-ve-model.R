###################################################################################################
#Title: Severe Waning VE Model
#Author: Hailey Park
#Date: March 25, 2023
###################################################################################################

rm(list=ls())

#Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)
library(scales)
library(lme4)
library(merTools)
library(sjstats)

#Load in raw datasets (these .csv files are from an excel doc of my manual review of recent waning literature)
waning_data <- read.csv("data/processed-data/severe waning data absolute 120324 mean.csv") 
weights <- read.csv("data/processed-data/severe waning data absolute 120324 weights.csv") 
waning_data_95ui <- read.csv("data/processed-data/severe waning data absolute 120324 95ui.csv")

#Reformat data to long and clean data
# NOTE: Mean and 95% CI data point estimates are taken directly from the literature. The weights are calculated directly 
#       from the 95% CI as (1/[upper 95% CI - lower 95% CI]), and used in the lmer model to weigh how much to use one 
#       study's estimates from another (e.g. VE estimates with narrower 95% CI are weighed greater than those with wider 95% CI)

long_data_mean <- melt(waning_data) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve =  value/100, 
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "Boosted"), ordered = TRUE),
         age_group = Age, 
         study = as.factor(References),
         estimate = "mean") %>%
  dplyr::select(c(months, prior_inf, ve, num_doses, age_group, study, estimate)) %>%
  mutate(ve_input = log(1 - ve),  #log of 1-VE
         month_input = log(months)) 

long_data_lower <- melt(waning_data_95ui %>% filter(X95UI == "lower")) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve =  value/100,
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "Boosted"), ordered = TRUE),
         age_group = Age, 
         study = as.factor(References),
         estimate = "lower") %>%
  dplyr::select(c(months, prior_inf, ve, num_doses, age_group, study, estimate)) %>%
  mutate(ve_input = log(1 - ve),  #log of 1-VE
         month_input = log(months))

long_data_upper <- melt(waning_data_95ui %>% filter(X95UI == "upper")) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve =  value/100, 
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "Boosted"), ordered = TRUE),
         age_group = Age, 
         study = as.factor(References),
         estimate = "upper") %>%
  dplyr::select(c(months, prior_inf, ve, num_doses, age_group, study, estimate)) %>%
  mutate(ve_input = log(1 - ve),  #log of 1-VE
         month_input = log(months))

long_data_weights <- melt(weights) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  #months
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "Boosted"), ordered = TRUE),
         age_group = Age, 
         study = as.factor(References),
         weight = value) %>%
  dplyr::select(c(months, prior_inf, num_doses, age_group, study, weight))


#Merge weights to raw waning data
merged_data <- merged_data <- merge(rbind(long_data_mean, long_data_lower, long_data_upper),
                                    long_data_weights, by = c("months", "prior_inf", "num_doses", "age_group", "study"), all.x = TRUE) %>%
  mutate(estimate = factor(estimate, levels = c("lower", "mean", "upper")))

  
#Fit model (The model is fit to specific immunity types -- infection-acquired, vaccine-derived, hybrid -- make sure to FILTER correct immunity type BEFOREHAND)
severe_model <- lmer(ve_input ~ month_input  + factor(age_group)  + prior_inf + estimate + (month_input|study),
                      weights = weight,
                      data = merged_data %>% filter(num_doses == 'Unvaccinated', prior_inf == 1)) #infection-acquired immunity
                      #data = merged_data %>% filter(num_doses == 'Boosted', prior_inf == 1)) #hybrid immunity
                      #data = merged_data %>% filter(num_doses == 'Boosted', prior_inf == 0)) #vaccine-derived immunity



summary(severe_model)


#Prediction for waning model
# NOTE: We are creating an empty dataframe to populate waning predictions. Since we are fitting the waning model by each immunity type,
#       we are creating predictions for each model/immunity type, so will need to select the correct "prior_inf" value according to immunity
#       type.
new_data_weekly <- data.frame(age_group = rep(c("18-49 years", "50-64 years", "65+ years"), each = 3),
                       estimate = rep(c("lower", "mean", "upper"), times = 3),
                       prior_inf = rep(c(1), each = 1, times = 9), #infection-acquired or hybrid immunity
                       #prior_inf = rep(c(0), each = 1, times = 9), #vaccine-derived immunity
                       study = "NA")
prediction_data_weekly <- new_data_weekly[rep(seq_len(nrow(new_data_weekly)), 104), ]
prediction_data_weekly$month_input <- rep(log(c(1:104)/4.345), each = 9)
prediction_data_weekly$months <- rep(floor(c(1:104)/4.345), each = 9)
prediction_data_weekly$weeks <- rep((c(1:104)), each = 9)

preds_weekly <- predict(severe_model, newdata = prediction_data_weekly, allow.new.levels = TRUE)
prediction_data_weekly$ve_pred <- 1 - exp(preds_weekly)

#Split the 65+ age group into a 65-74 and 75+ age group
below65 <- prediction_data_weekly %>% filter(age_group != "65+ years")
above65_74 <- prediction_data_weekly %>% filter(age_group == "65+ years") %>%
  mutate(age_group = "65-74 years")
above75 <- prediction_data_weekly %>% filter(age_group == "65+ years") %>%
  mutate(age_group = "75+ years")

#Combine all age groups together, and filter out negative VEs
prediction_data_immunocompetent <- rbind(below65, above65_74, above75) %>%
  rowwise() %>% mutate(ve_pred = max(0, ve_pred))

#Write immunocompetent curves to csv (make sure to select correct immunity type file to write)
write.csv(prediction_data_immunocompetent, "data/clean-data/severe_waning_predictions_weekly_infectionImmunityOnly.csv")
# write.csv(prediction_data_immunocompetent, "data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly.csv")
# write.csv(prediction_data_immunocompetent, "data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly.csv")

#Create immunocompromised curves (15% decrement) and save to csv
# NOTE: For infection-acquired and vaccine-derived immunity, doing a 15% overall decrement to immunocompetent curve.
#       For hybrid immunity, I did a 5% decrement to the immunocompetent curve, plus additional 10% decrement through adjustment to the waning rate.
#       This was done to better match the relative ve literature estimates for hybrid immunity in immunocompromised populations

infection_immunocompromised <- read.csv("data/clean-data/severe_waning_predictions_weekly_infectionImmunityOnly.csv")[,-1] %>% rowwise() %>%
  mutate(ve_pred = max(0, ve_pred - 0.15))
vaccine_immunocompromised <- read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1] %>% rowwise() %>%
  mutate(ve_pred = max(0, ve_pred - 0.15))
hybrid_immunocompromised <- read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1] %>% rowwise() %>%
  mutate(ve_pred = max(0, ve_pred - (weeks * (0.00416/4.345 * 2)) - 0.05))
  

# write.csv(infection_immunocompromised, "data/clean-data/severe_waning_predictions_weekly_infectionImmunityOnly_immunocompromised.csv")
# write.csv(vaccine_immunocompromised, "data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly_immunocompromised.csv")
# write.csv(hybrid_immunocompromised, "data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")

#Plot curves (read in correct waning curve file based on immunity type and immunocompromised status)
plot_data <- read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")[,-1]


plot_data %>% 
  filter(age_group == "65-74 years") %>%
  ggplot(aes(weeks, ve_pred, group=estimate, color=estimate)) + 
  geom_line() + geom_point() +
  ylim(0,1) +
  ylab("Protective Effectiveness (%)") +
  xlab("Time Since (months)") +
  ggtitle("Hybrid immunity waning curves against severe COVID-19\nAge Group: 65+ years\nRisk Group: Immunocompromised")


