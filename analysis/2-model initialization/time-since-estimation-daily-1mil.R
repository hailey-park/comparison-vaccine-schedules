###################################################################################################
#Title: 'Time Since Last' Estimation
#Author: Hailey Park
#Date: June 19, 2024
###################################################################################################

rm(list=ls())

#setwd(here::here())

#Load libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)
library(doParallel)
library(foreach)

#Detect cores
num_cores <- detectCores()
registerDoParallel(num_cores - 1)

#Read data
cases_by_day <- read.csv("data/clean-data/cases_by_day.csv")[,-1]
booster_doses_1st_by_day <- read.csv("data/clean-data/booster_1st_doses_by_day.csv")[,-1]
booster_doses_2nd_by_day <- read.csv("data/clean-data/booster_2nd_doses_by_day_updated.csv")[,-1]
fully_vax_doses_by_day <- read.csv("data/clean-data/fully_vax_doses_by_day.csv")[,-1]

#'Time Since Last' Estimation -- Initializing model at December 31, 2021
#NOTE: Estimation run separately for each age group because of differences in vaccine coverage and seroprevalence

#Clean data
booster_doses_1st_by_day$day <- as.character(booster_doses_1st_by_day$day)
booster_doses_2nd_by_day$day <- as.character(booster_doses_2nd_by_day$day)
fully_vax_doses_by_day$day <- as.character(fully_vax_doses_by_day$day)
cases_by_day$day <- as.character(cases_by_day$date)

#Calculate time since last vaccine dose or infection
add.days= function(date,n) {seq(date, by = paste (n, "days"), length = 2)[2]}

time_since_last <- function(df) {
  
  #Get age-specific dose/case distributions
  booster_doses_1st_by_day <- booster_doses_1st_by_day %>% filter(age_group == df$age_group[1])
  booster_doses_2nd_by_day <- booster_doses_2nd_by_day %>% filter(age_group == df$age_group[1])
  fully_vax_doses_by_day <- fully_vax_doses_by_day %>% filter(age_group == df$age_group[1])

  #Calculate time since last dose and time since last infection
  set.seed(88)
  last_dose_and_inf <- df %>% mutate(time_since_last_dose = ifelse(prior_vacc == 'fullvax', sample(as.character(fully_vax_doses_by_day$day), size = sum(prior_vacc == 'fullvax'), prob = fully_vax_doses_by_day$perc_doses, replace = TRUE),
                                                                   ifelse(prior_vacc == 'boosted_1st', sample(as.character(booster_doses_1st_by_day$day), size = sum(prior_vacc == 'boosted_1st'), prob = booster_doses_1st_by_day$perc_doses, replace = TRUE),
                                                                          ifelse(prior_vacc == 'boosted_2nd', sample(as.character(booster_doses_2nd_by_day$day), size = sum(prior_vacc == 'boosted_2nd'), prob = booster_doses_2nd_by_day$perc_doses, replace = TRUE),
                                                                                 NA))),
                                     time_since_last_inf = ifelse(prior_inf %in% c('inf', 'reinf'), sample(as.character(cases_by_day$day), size = sum(prior_inf %in% c('inf', 'reinf')), prob = cases_by_day$perc_cases, replace = TRUE),
                                                                  NA)) 
  
  reinf_only <- last_dose_and_inf %>% filter(prior_inf == 'reinf')
  
#system.time({
  # Return a data frame
  reinf_estimation <- rbind(foreach (i=1:(floor(nrow(reinf_only)/10000) - 1), .combine=rbind) %dopar% {
    
    reinf_only[(1 + (10000 * (i-1))):(10000 * i),] %>%
        rowwise() %>% mutate(time_since_last_reinf = sample(as.character(cases_by_day$day[as.Date(cases_by_day$day) >= as.Date(time_since_last_inf)]),
                                                            size = 1,
                                                            prob = cases_by_day$perc_cases[as.Date(cases_by_day$day) >= as.Date(time_since_last_inf)],
                                                            replace = TRUE))
  },
  reinf_only[(1 + (10000 * (floor(nrow(reinf_only)/10000) - 1))):nrow(reinf_only),] %>%
    rowwise() %>% mutate(time_since_last_reinf = sample(as.character(cases_by_day$day[as.Date(cases_by_day$day) >= as.Date(time_since_last_inf)]),
                                                        size = 1,
                                                        prob = cases_by_day$perc_cases[as.Date(cases_by_day$day) >= as.Date(time_since_last_inf)],
                                                        replace = TRUE)))

  combined_time_since <- merge(last_dose_and_inf, reinf_estimation, all.x = TRUE) %>%
    mutate(time_since_last_dose_inf = pmax(as.Date(time_since_last_dose), as.Date(time_since_last_inf), as.Date(time_since_last_reinf), na.rm =  TRUE))
  
  return(combined_time_since)
}

#Create matrices for each age group

set.seed(88)
age_0_17_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                           prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                           ntimes=c(0.08601, 0.332, 0.121, 0.30001, 0.066 * (1-0.476), 0.066 * 0.47601, 0.074 * (1-0.476), 0.074 * 0.47601, 0.009 * (1-0.476), 0.009 * 0.47601, 0.012) * (1000000*0.2176))

age_0_17_cal <- as.data.frame(lapply(age_0_17_specs[,1:2], rep, age_0_17_specs$ntimes)) %>%
  mutate(individual = c(1:(1000000*0.2176)),
         age_group = "0-17 years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.2176), prob = c(0.882, 0.026, 0.092), replace = TRUE))
  
 
set.seed(88)
age_18_29_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                             prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                             ntimes=c(0.04601, 0.17701, 0.1, 0.246, 0.175 * (1-0.388), 0.175 * 0.388, 0.194 * (1-0.388), 0.194 * 0.38801, 0.043 * (1-0.388), 0.043 * 0.388, 0.019) * (1000000*0.157))

age_18_29_cal <- as.data.frame(lapply(age_18_29_specs[,1:2], rep, age_18_29_specs$ntimes)) %>%
  mutate(individual = c((1 + (1000000*0.2176)):((1000000*0.3746))),
         age_group = "18-29 years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.157), prob = c(0.726, 0.033, 0.241), replace = TRUE))

set.seed(88)
age_30_49_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                              prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                              ntimes=c(0.02601, 0.10101, 0.074, 0.185, 0.232 * (1-0.271), 0.232 * 0.271, 0.258 * (1-0.271), 0.258 * 0.271005, 0.1 * (1-0.271), 0.1 * 0.271, 0.024) * (1000000*0.2621))

age_30_49_cal <- as.data.frame(lapply(age_30_49_specs[,1:2], rep, age_30_49_specs$ntimes)) %>%
  mutate(individual = c((1 + (1000000*0.3746)):((1000000*0.6367))),
         age_group = "30-49 years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.2621), prob = c(0.726, 0.056, 0.218), replace = TRUE))

set.seed(88)
age_50_64_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                              prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                              ntimes=c(0.02201, 0.08501, 0.048005, 0.118, 0.274 * (1-0.420), 0.274 * 0.420006, 0.304 * (1-0.420), 0.304 * 0.420, 0.131 * (1-0.420), 0.131 * 0.420, 0.018) * (1000000*0.1868))

age_50_64_cal <- as.data.frame(lapply(age_50_64_specs[,1:2], rep, age_50_64_specs$ntimes)) %>%
  mutate(individual = c((1 + (1000000*0.6367)):((1000000*0.8235))),
         age_group = "50-64 years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.1868), prob = c(0.366, 0.089, 0.545), replace = TRUE))

set.seed(88)
age_65_74_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                              prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                              ntimes=c(0.01601, 0.06002, 0.035, 0.086, 0.272 * (1-0.686), 0.272 * 0.68601, 0.303 * (1-0.686), 0.303 * 0.686, 0.208 * (1-0.686), 0.208 * 0.686, 0.02) * (1000000*0.1035))

age_65_74_cal <- as.data.frame(lapply(age_65_74_specs[,1:2], rep, age_65_74_specs$ntimes)) %>%
  mutate(individual = c((1 + (1000000*0.8235)):((1000000*0.927))),
         age_group = "65-74 years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.1035), prob = c(0.124, 0.1065, 0.7695), replace = TRUE))

set.seed(88)
age_75plus_specs <- data.frame(prior_vacc=c("unvax","unvax","fullvax", "fullvax", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "boosted_1st", "boosted_2nd", "fullvax"), 
                              prior_inf=c("inf","reinf","inf","reinf","inf", "inf", "reinf", "reinf", "noinf", "noinf", "noinf"), 
                              ntimes=c(0.01601, 0.06002, 0.03501, 0.08603, 0.272 * (1-0.736), 0.272 * 0.73601, 0.303 * (1-0.736), 0.303 * 0.736008, 0.208 * (1-0.736), 0.208 * 0.736, 0.02) * (1000000*0.073))

age_75plus_cal <- as.data.frame(lapply(age_75plus_specs[,1:2], rep, age_75plus_specs$ntimes)) %>%
  mutate(individual = c((1 + (1000000*0.927)):(1000000)),
         age_group = "75+ years",
         risk_group = sample(c('healthy', 'immunocompromised', 'higher risk'), (1000000*0.073), prob = c(0.124, 0.1065, 0.7695), replace = TRUE))


set.seed(88)
time_since_results <- list(age_0_17_cal, age_18_29_cal, age_30_49_cal, age_50_64_cal, age_65_74_cal, age_75plus_cal) %>%
  lapply(time_since_last)

inspection <- time_since_results[[5]] 



#combine the age groups into one dataframe
combined <- bind_rows(time_since_results) 


write.csv(combined, "data/clean-data/entire_population_model_initialization_daily_1mil_updated.csv")
############################################################
#PLOTS

pdf("time-since-plots-1mil.pdf")

inspection <- read.csv("data/clean-data/entire_population_model_initialization_daily_1mil_updated.csv")[,-1] %>%
  filter(age_group == "0-17 years") %>% dplyr::select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
  group_by(prior_vacc, time_since_last_dose) %>% summarise(total = n()) %>% filter(prior_vacc != 'unvax')

for(i in c(1:6)) {
  age_group_label <- case_when(i == 1 ~ "0-17 years",
                               i == 2 ~ "18-29 years",
                               i == 3 ~ "30-49 years",
                               i == 4 ~ "50-64 years",
                               i == 5 ~ "65-74 years",
                               i == 6 ~ "75+ years",
                               TRUE ~ NA) 
  
  inspection <- time_since_results[[i]] %>% dplyr::select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
    group_by(prior_vacc, time_since_last_dose) %>% summarise(total = n()) %>% filter(prior_vacc != 'unvax')
  
  observed_data <- as.data.frame(rbind(fully_vax_doses_by_day %>% mutate(prior_vacc = "observed_fullvax"),
                                       booster_doses_1st_by_day %>% mutate(prior_vacc = "observed_boosted_1st"),
                                       booster_doses_2nd_by_day %>% mutate(prior_vacc = "observed_boosted_2nd"))) %>% 
    dplyr::select(prior_vacc, day, age_group, perc_doses) %>% filter(age_group == age_group_label)
  
  #plot simulated prior vacc
  plot(ggplot(data = inspection, aes(x = as.Date(time_since_last_dose), y = total, color = prior_vacc)) +
         geom_line() +
         ylab("Total COVID-19 Doses") +
         xlab("days") +
         scale_x_date(date_labels = ("%b-%Y"),
                      limits = c(as.Date("2020-12-15"), as.Date("2023-07-01")),
                      breaks = "1 months") +
         ylim(0, NA)+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         ggtitle(paste0("Simulated Time-Since Last Vaccination\nAge Group: ", age_group_label)))
  
  
  inspection <- time_since_results[[i]] %>% dplyr::select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
    group_by(time_since_last_inf) %>% summarise(total = n()) #%>% filter(!is.na(time_since_last_reinf))
  
  observed_data <- cases_by_day 
  
  ## plot simulated prior inf
  plot(ggplot(data = inspection, aes(x = as.Date(time_since_last_inf), y = total)) +
         geom_line() +
         ylab("Total COVID-19 Infections") +
         xlab("Time") +
         scale_x_date(date_labels = ("%b-%Y"),
                      limits = c(as.Date("2020-03-01"), as.Date("2023-07-01")),
                      breaks = "1 months") +
         ylim(0, NA)+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         ggtitle(paste0("Simulated Time-Since Last Infection\nAge Group: ", age_group_label)))
  
  inspection <- time_since_results[[i]] %>% dplyr::select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
    group_by(time_since_last_inf) %>% summarise(total = n()) %>% filter(!is.na(time_since_last_reinf))
  
  observed_data <- cases_by_day 
  
  ## plot simulated prior reinf
  plot(ggplot(data = inspection, aes(x = as.Date(time_since_last_reinf), y = total)) +
         geom_line() +
         ylab("Total COVID-19 Infections") +
         xlab("Time") +
         scale_x_date(date_labels = ("%b-%Y"),
                      limits = c(as.Date("2020-03-01"), as.Date("2023-07-01")),
                      breaks = "1 months") +
         ylim(0, NA)+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         ggtitle(paste0("Simulated Time-Since Last Re-infection\nAge Group: ", age_group_label)))
  
  
  inspection <- time_since_results[[i]] %>% dplyr::select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_dose_inf) %>%
    group_by(time_since_last_dose_inf) %>% summarise(total = n()) %>%
    mutate(time_since_adjusted = floor_date(time_since_last_dose_inf, unit = "day")) %>%
    group_by(time_since_adjusted) %>% summarise(total = sum(total))
  
  
  #plot simulated time since
  plot(ggplot(data = inspection, aes(x = as.Date(time_since_last_dose_inf), y = total)) +
         geom_line() +
         ylab("Total COVID-19 Immune Events") +
         xlab("Time") +
         scale_x_date(date_labels = ("%b-%Y"),
                      limits = c(as.Date("2021-01-01"), as.Date("2023-07-01")),
                      breaks = "1 months") +
         ylim(0, NA)+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         ggtitle(paste0("Simulated Time-Since Last Immune Event\nAge Group: ", age_group_label)))
  
}


dev.off()