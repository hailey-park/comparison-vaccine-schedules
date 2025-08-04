########################################################################################################################
#Title: Tables for primary outcomes (vaccination strategy comparison)
#Author: Hailey Park
#Date: April 20, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(reshape2)
library(data.table)

#cleaning results (general risk group-specific)
cleaning_sim_general <- function(df) {
  clean_df <- df %>%
    mutate(severe_total = select(., day1:day547) %>% rowSums(na.rm = TRUE),
           nonsevere_total = select(., nonsevere_day1:nonsevere_day547) %>% rowSums(na.rm = TRUE)) %>%
    dplyr::select(age_group, risk_group, total_pop, severe_total, nonsevere_total, total_vaccines, total_vaccinated, total_severe_vax, total_nonsevere_vax) %>%
    dplyr::mutate(general_risk_cat = case_when(age_group %in% c("65-74 years", "75+ years") ~ "65+ years",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "healthy" ~ "18-64 years, healthy",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "higher risk" ~ "18-64 years, higher risk",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "immunocompromised" ~ "18-64 years, immunocompromised",
                                        TRUE ~ age_group)) %>%
    dplyr::group_by(general_risk_cat) %>% dplyr::summarise(across(total_pop:total_nonsevere_vax, sum)) %>%
    dplyr::mutate(severe_risk_vax = (total_severe_vax/total_vaccinated)/1.5 * 100000)
  
  return(clean_df)
  
}


#Read in strategies
setwd("~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/")
strat_real <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/'), pattern = "upper", invert = FALSE, value = TRUE), pattern = ".csv", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)

setwd("~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic")
strat_0 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_0", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)

setwd("~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic")
strat_1 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_1", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_2 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_2", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_3 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_3", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_4 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_4", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_5 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_5", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_6 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_6", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_7 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_7", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_8 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_8", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_9 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_9", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)
strat_10 <- lapply(grep(grep(list.files('~/Stanford Research/acip-analysis/sherlock-testing/final-analysis/simulation-results-11/realistic'), pattern = "lower|upper", invert = TRUE, value = TRUE), pattern = "strat_10", value = TRUE), read.csv, header=T, sep=',', dec='.')%>%
  lapply(cleaning_sim_general)



ui_95_lower  <- function(x) {
  round(mean(x) + (-1.96) * sd(x)/sqrt(length(x)), 1)
}


ui_95_upper  <- function(x) {
  round(mean(x) + (1.96) * sd(x)/sqrt(length(x)), 1)
}


cleaned_results_mean <- list(do.call(rbind, strat_real) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        do.call(rbind, strat_0) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        do.call(rbind, strat_1) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        do.call(rbind, strat_2) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        #do.call(rbind, strat_3) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        do.call(rbind, strat_4) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        do.call(rbind, strat_6) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                         do.call(rbind, strat_5) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        #  do.call(rbind, strat_7) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) , 
                        # do.call(rbind, strat_8) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE) ,
                      do.call(rbind, strat_9) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE))
                       #do.call(rbind, strat_10) %>% group_by(general_risk_cat) %>% summarise_all(mean, na.rm = TRUE))

cleaned_results_lower <- list(do.call(rbind, strat_real) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_0) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_1) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_2) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             #do.call(rbind, strat_3) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_4) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_6) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_5) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             # do.call(rbind, strat_7) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             # do.call(rbind, strat_8) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) , 
                             do.call(rbind, strat_9) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower) ) 
                             #do.call(rbind, strat_10) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_lower)) 


cleaned_results_upper <- list(do.call(rbind, strat_real) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_0) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_1) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_2) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             #do.call(rbind, strat_3) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_4) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_6) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_5) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             # do.call(rbind, strat_7) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             # do.call(rbind, strat_8) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) , 
                             do.call(rbind, strat_9) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper) ) 
                             #do.call(rbind, strat_10) %>% group_by(general_risk_cat) %>% summarise_all(ui_95_upper)) 



processing_sim_mean <- function(df) {
  
  reference_group <- cleaned_results_mean[[2]] %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  processed_df <- df %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  final_df <- processed_df %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, total_vaccines, total_vaccinated, total_severe_vax, total_pop) %>%
    ungroup() %>%
    mutate(ARR = reference_group$severe_risk - severe_risk_vax,
           RRR = (reference_group$severe_risk - severe_risk)/reference_group$severe_risk,
           NNT = 1/((ARR)/100000)) %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, ARR, NNT, RRR, total_vaccines)
  
  return(final_df)
}

processing_sim_lower <- function(df) {
  
  reference_group <- cleaned_results_lower[[2]] %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  processed_df <- df %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  final_df <- processed_df %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, total_vaccines, total_vaccinated, total_severe_vax) %>%
    ungroup() %>%
    mutate(ARR = reference_group$severe_risk - severe_risk_vax,
           RRR = (reference_group$severe_risk - severe_risk)/reference_group$severe_risk,
           NNT = 1/((ARR)/100000)) %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, ARR, NNT, RRR, total_vaccines)
  
  return(final_df)
}


processing_sim_upper <- function(df) {
  
  reference_group <- cleaned_results_upper[[2]] %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  processed_df <- df %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000)
  
  final_df <- processed_df %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, total_vaccines, total_vaccinated, total_severe_vax) %>%
    ungroup() %>%
    mutate(ARR = reference_group$severe_risk - severe_risk_vax,
           RRR = (reference_group$severe_risk - severe_risk)/reference_group$severe_risk,
           NNT = 1/((ARR)/100000)) %>%
    dplyr::select(general_risk_cat, severe_total, severe_risk, severe_risk_vax, ARR, NNT, RRR, total_vaccines)
  
  return(final_df)
}



cleaned_results_mean[c(1:8)] %>% lapply(processing_sim_mean)
cleaned_results_lower[c(1:8)] %>% lapply(processing_sim_lower)
cleaned_results_upper[c(1:8)] %>% lapply(processing_sim_upper)


