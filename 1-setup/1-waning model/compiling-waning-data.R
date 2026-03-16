###################################################################################################
#Title: Compiling waning curves for model input
#Author: Hailey Park
#Date: July 30, 2025
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


#Read in all waning curves: by immunity type (infection-acquired, vaccine-derived, hybrid), severity (severe, general infection), and immunocompromised status
severe_waning_inf <- read.csv("data/clean-data/severe_waning_predictions_weekly_infectionImmunityOnly.csv")[,-1]
severe_waning_inf_immuno <- read.csv("data/clean-data/severe_waning_predictions_weekly_infectionImmunityOnly_immunocompromised.csv")[,-1]

severe_waning_vacc <- read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]
severe_waning_vacc_immuno <- read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly_immunocompromised.csv")[,-1]

severe_waning_hybrid <- read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]
severe_waning_hybrid_immuno <- read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")[,-1]

nonsevere_waning_inf <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_infectionImmunityOnly.csv")[,-1]
nonsevere_waning_inf_immuno <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_infectionImmunityOnly_immunocompromised.csv")[,-1]

nonsevere_waning_vacc <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]
nonsevere_waning_vacc_immuno <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_vaccineImmunityOnly_immunocompromised.csv")[,-1]

nonsevere_waning_hybrid <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]
nonsevere_waning_hybrid_immuno <- read.csv("data/clean-data/nonsevere_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")[,-1]

#Cleaning function (make sure to select correct estimate -- mean, or lower/upper 95% CI estimates)
df_cleaning_func <- function(x) { 
  x %>%
    filter(estimate == 'upper') %>% ##SELECT CORRECT ESTIMATE HERE (mean, lower, upper)
    dplyr::select(-any_of(c("study", "estimate", 'month_input', "months", "prior_inf"))) %>%
    rowwise() %>% mutate(ve_pred = max(ve_pred, 0)) %>% 
    dplyr::group_by(weeks) %>% dplyr::summarise(ve_pred = mean(ve_pred))
}



#Compile all waning curves into 1 dataframe. Make sure to select correct estimate (mean, lower, upper) in the cleaning function and when you write file to csv.
# NOTE: I am writing 2 waning curve files, one is this section, the second is the "alt" file. This first section creates a wide waning curve data file. 
#       The second creating a long waning curve file. I use both in the model inputs.

#First waning curve df (wide format)
waning_curves <- list(severe_waning_inf, severe_waning_inf_immuno, severe_waning_vacc, severe_waning_vacc_immuno, severe_waning_hybrid, severe_waning_hybrid_immuno, nonsevere_waning_inf, nonsevere_waning_inf_immuno, nonsevere_waning_vacc, nonsevere_waning_vacc_immuno, nonsevere_waning_hybrid, nonsevere_waning_hybrid_immuno) %>%
  lapply(df_cleaning_func)

severe_waning_inf <- waning_curves[[1]] %>% mutate(immuno = 0) %>% dplyr::rename(severe_inf_ve = ve_pred) %>% arrange(weeks)
severe_waning_inf_immuno <- waning_curves[[2]] %>% mutate(immuno = 1) %>% dplyr::rename(severe_inf_ve = ve_pred) %>% arrange(weeks)
severe_waning_vacc <- waning_curves[[3]] %>% mutate(immuno = 0) %>% dplyr::rename(severe_vacc_ve = ve_pred) %>% arrange(weeks)
severe_waning_vacc_immuno <- waning_curves[[4]] %>% mutate(immuno = 1) %>% dplyr::rename(severe_vacc_ve = ve_pred) %>% arrange(weeks)
severe_waning_hybrid <- waning_curves[[5]] %>% mutate(immuno = 0) %>% dplyr::rename(severe_hybrid_ve = ve_pred) %>% arrange(weeks)
severe_waning_hybrid_immuno <- waning_curves[[6]] %>% mutate(immuno = 1) %>% dplyr::rename(severe_hybrid_ve = ve_pred) %>% arrange(weeks)

nonsevere_waning_inf <- waning_curves[[7]] %>% mutate(immuno = 0) %>% dplyr::rename(nonsevere_inf_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_inf_immuno <- waning_curves[[8]] %>% mutate(immuno = 1) %>% dplyr::rename(nonsevere_inf_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_vacc <- waning_curves[[9]] %>% mutate(immuno = 0) %>% dplyr::rename(nonsevere_vacc_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_vacc_immuno <- waning_curves[[10]] %>% mutate(immuno = 1) %>% dplyr::rename(nonsevere_vacc_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_hybrid <- waning_curves[[11]] %>% mutate(immuno = 0) %>% dplyr::rename(nonsevere_hybrid_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_hybrid_immuno <- waning_curves[[12]] %>% mutate(immuno = 1) %>% dplyr::rename(nonsevere_hybrid_ve = ve_pred) %>% arrange(weeks)


waning_data_clean <- merge(merge(merge(merge(merge(rbind(severe_waning_inf, severe_waning_inf_immuno), rbind(severe_waning_vacc, severe_waning_vacc_immuno), by = c("immuno", "weeks")),
                                             rbind(severe_waning_hybrid, severe_waning_hybrid_immuno), by = c("weeks", "immuno"), all.x = TRUE),
                                       rbind(nonsevere_waning_inf, nonsevere_waning_inf_immuno), by = c("weeks", "immuno"), all.x = TRUE),
                                 rbind(nonsevere_waning_vacc, nonsevere_waning_vacc_immuno), by = c("weeks", "immuno"), all.x = TRUE),
                           rbind(nonsevere_waning_hybrid, nonsevere_waning_hybrid_immuno), by = c("weeks", "immuno"), all.x = TRUE)


write.csv(waning_data_clean, "data/clean-data/waning_data_clean_upper.csv")
write.csv(waning_data_clean, "data/clean-data/waning_data_clean_lower.csv")
write.csv(waning_data_clean, "data/clean-data/waning_data_clean.csv") #MEAN


#Second waning curve df (long format)
severe_waning_inf <- waning_curves[[1]] %>% mutate(immuno = 0, prior_inf = 1, prior_vacc = 0) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)
severe_waning_inf_immuno <- waning_curves[[2]] %>% mutate(immuno = 1, prior_inf = 1, prior_vacc = 0) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)
severe_waning_vacc <- waning_curves[[3]] %>% mutate(immuno = 0, prior_inf = 0, prior_vacc = 1) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)
severe_waning_vacc_immuno <- waning_curves[[4]] %>% mutate(immuno = 1, prior_inf = 0, prior_vacc = 1) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)
severe_waning_hybrid <- waning_curves[[5]] %>% mutate(immuno = 0, prior_inf = 1, prior_vacc = 1) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)
severe_waning_hybrid_immuno <- waning_curves[[6]] %>% mutate(immuno = 1, prior_inf = 1, prior_vacc = 1) %>% dplyr::rename(severe_ve = ve_pred) %>% arrange(weeks)

nonsevere_waning_inf <- waning_curves[[7]] %>% mutate(immuno = 0, prior_inf = 1, prior_vacc = 0) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_inf_immuno <- waning_curves[[8]] %>% mutate(immuno = 1, prior_inf = 1, prior_vacc = 0) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_vacc <- waning_curves[[9]] %>% mutate(immuno = 0, prior_inf = 0, prior_vacc = 1) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_vacc_immuno <- waning_curves[[10]] %>% mutate(immuno = 1, prior_inf = 0, prior_vacc = 1) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_hybrid <- waning_curves[[11]] %>% mutate(immuno = 0, prior_inf = 1, prior_vacc = 1) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)
nonsevere_waning_hybrid_immuno <- waning_curves[[12]] %>% mutate(immuno = 1, prior_inf = 1, prior_vacc = 1) %>% dplyr::rename(nonsevere_ve = ve_pred) %>% arrange(weeks)


waning_data_clean <- merge(rbind(severe_waning_inf, severe_waning_inf_immuno, severe_waning_vacc, severe_waning_vacc_immuno, severe_waning_hybrid, severe_waning_hybrid_immuno),
                           rbind(nonsevere_waning_inf, nonsevere_waning_inf_immuno, nonsevere_waning_vacc, nonsevere_waning_vacc_immuno, nonsevere_waning_hybrid, nonsevere_waning_hybrid_immuno),
                           by = c("weeks", "immuno", "prior_inf", "prior_vacc"), all.x = TRUE)


write.csv(waning_data_clean, "data/clean-data/waning_data_clean_alt_upper.csv")
write.csv(waning_data_clean, "data/clean-data/waning_data_clean_alt_lower.csv")
write.csv(waning_data_clean, "data/clean-data/waning_data_clean_alt.csv") #MEAN

