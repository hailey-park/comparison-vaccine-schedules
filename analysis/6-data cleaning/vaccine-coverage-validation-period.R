###################################################################################################
#Title: Cleaning Vaccine Coverage Data Used for Model Validation (July 1, 2024 - Jan 1, 2025)
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

#Read in data
  #Vaccine data from  October 2024 - March 2025 from this website (https://data.cdc.gov/Vaccinations/Weekly-Cumulative-COVID-19-Vaccination-Coverage-an/ksfb-ug5d/about_data)
  #Vaccine data for children from  October 2024 - March 2025 from this website (https://data.cdc.gov/Child-Vaccinations/Weekly-Parental-Intent-for-Vaccination-and-Cumulat/ker6-gs6z/about_data)
vaccine_2025 <- read.csv("raw-data/Weekly_Cumulative_COVID-19_Vaccination_Coverage_and_Intent__Overall__by_Selected_Demographics_and_Jurisdiction__Among_Adults_18_Years_and_Older_20250321.csv")
vaccine_2025_child <- read.csv("raw-data/Weekly_Parental_Intent_for_Vaccination_and_Cumulative_Percentage_of_Children_6_Months_-17_Years_Who_are_Up_to_date_with_the_COVID-19_Vaccines_by_Season__United_States_20250321.csv")

  #Vaccine data from September 2024- March 2025 from this website (https://data.cdc.gov/Vaccinations/National-Immunization-Survey-Adult-COVID-Module-NI/si7g-c2bs/about_data)
vaccines_risk <- read.csv("raw-data/National_Immunization_Survey_Adult_COVID_Module__NIS-ACM___RespVaxView__Data___Centers_for_Disease_Control_and_Prevention__cdc.gov__20250321.csv")

  #US population prevalence estimates by risk group and age group
prevalence <- read.csv("data/processed-data/risk_group_prevalence_estimates_1+comorbidity.csv")


#Clean data

  #Prevalence with adjusted risk groups
prevalence_adj <- prevalence %>% mutate(high_risk = if_else(risk_group != "healthy", "Yes", "No")) %>%
  group_by(age_group, high_risk) %>% summarise(prevalence = sum(prevalence)) %>%
  group_by(age_group) %>% mutate(prop_prev_by_age = prevalence/(sum(prevalence))) %>%
  mutate(age_group = case_when(age_group == "0-17 years (Children)" ~ "0-17 years",
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) %>%
  dplyr::select(age_group, high_risk, prop_prev_by_age) %>%
  pivot_wider(names_from = high_risk,
              values_from = prop_prev_by_age)
  
  
  #Vaccine uptake by risk (healthy vs. composite high risk (includes our defined 'higher risk' and 'immunocompromised' groups))
vaccines_by_risk <- vaccines_risk %>% filter(Geography.Type == "National",
                                             Group.Name == "Health condition associated with higher risk for viral respiratory diseases",
                                             Indicator.Name == "Vaccination coverage and intent (among all adults 18+)",
                                             Indicator.Category == "Received a 2024-2025 COVID-19 dose ") %>%
  mutate(month = case_when(Time.Period == "September 1 - September 28" ~ as.Date("2024-09-01"),
                           Time.Period == "September 29 - October 26 " ~ as.Date('2024-10-01'),
                           Time.Period == "October 27 - November 30" ~ as.Date('2024-11-01'),
                           Time.Period == "December 1 - December 28 " ~ as.Date('2024-12-01'),
                           Time.Period == "December 29 - January 25" ~ as.Date("2025-01-01"),
                           Time.Period == "January 26 - February 22" ~ as.Date('2025-02-01'),
                           TRUE ~ NA)) %>%
  dplyr::select(c(Group.Category, month, Estimate....)) %>% rename(high_risk = Group.Category, vaccine_coverage = Estimate....) %>%
  unique() %>%
  pivot_wider(names_from = high_risk,
              values_from = vaccine_coverage) %>%
  mutate(ratio = Yes/No)

average_ratio <- mean(vaccines_by_risk$ratio)

#Vaccine doses by month from October 2024 through January 2025 (ADULTS)
doses_by_month_mar2025_adult <- vaccine_2025 %>% filter(Geographic.Level == "National", Demographic.Level == "Age", indicator_label == "Up-to-date", Demographic.Name %in% c("18-29 years", "30-39 years", "40-49 years", "50-64 years", "65-74 years", "75+ years")) %>%
  mutate(date = as.Date(Week_ending, format = "%m/%d/%Y"),
         age_group = if_else(Demographic.Name %in% c("30-39 years", "40-49 years"), "30-49 years", Demographic.Name)) %>% group_by(date, age_group) %>%
  summarise(perc = mean(Estimate)) %>%
  group_by(age_group) %>% mutate(perc_up_to_date_week = c((perc - lag(perc, 1))[2], (perc - lag(perc, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, perc, perc_up_to_date_week)) %>%
  filter(date >= as.Date("2024-06-29") & date <= as.Date("2025-01-01")) %>% dplyr::select(date, age_group, perc_up_to_date_week) %>%
  group_by(age_group)
  
#Vaccine doses by month from October 2024 through January 2025 (CHILDREN)
doses_by_month_mar2025_child <- vaccine_2025_child %>% filter(Geographic.Level == "National", Demographic_Level == "Age", Indicator_label == "Up-to-date", Demographic.Name %in% c("6 months-4 years", "5-17 years")) %>%
  mutate(date = as.Date(Week_ending)) %>% group_by(date) %>%
  summarise(perc = weighted.mean(Estimate, Unweighted.Sample.Size)) %>%
  mutate(perc_up_to_date_week = c((perc - lag(perc, 1))[2], (perc - lag(perc, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, perc, perc_up_to_date_week),
         age_group = "0-17 years") %>%
  filter(date >= as.Date("2024-7-01") & date <= as.Date("2025-01-01")) %>% dplyr::select(date, age_group, perc_up_to_date_week) %>%
  group_by(age_group) 

#Impute missing children's coverage data with gradual increase to match uptake in October
missing_children_coverage <- data.frame(date = (seq(as.Date('2024-06-30'), as.Date("2024-09-29"), by = "weeks")),
                                   age_group = rep(c("0-17 years"), each = 14),
                                   perc_up_to_date_week = 3.3/14) 

#Impute missing adult coverage data in September
missing_adult_coverage <- data.frame(date = rep(as.Date(c("2024-08-24", "2024-08-31")), 5),
                                     age_group = rep(c("18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"), each = 2),
                                     perc_up_to_date_week = c(rep((1.3 - 0.6)/3, 2), rep((1.55 - 0.6)/3, 2), rep((2.6 - 0.60755)/3, 2), rep((5.8 - 0.8)/3, 2), rep((4.6 - 1.6)/3, 2))) 

#Combine vaccine uptake data
combined <- rbind(missing_children_coverage, missing_adult_coverage, doses_by_month_mar2025_adult, doses_by_month_mar2025_child) %>%
  group_by(age_group) %>%
  arrange(date) %>%
  mutate(coverage = cumsum(perc_up_to_date_week))

#Calculate healthy group vaccine uptake
healthy_vax_coverage <- merge(combined %>% mutate(risk_group = "healthy"),
                              prevalence_adj, by = "age_group", all.x = TRUE) %>%
  mutate(risk_coverage = coverage/(Yes * average_ratio + No),
         risk_coverage_per_week = perc_up_to_date_week/(Yes * average_ratio + No))

#Add risk group information
combined_with_risk <- rbind(healthy_vax_coverage,
                            healthy_vax_coverage %>% mutate(risk_group = "immunocompromised",
                                                            risk_coverage = risk_coverage * average_ratio,
                                                            risk_coverage_per_week = risk_coverage_per_week * average_ratio),
                            healthy_vax_coverage %>% mutate(risk_group = "higher risk",
                                                            risk_coverage = risk_coverage * average_ratio,
                                                            risk_coverage_per_week = risk_coverage_per_week * average_ratio)) %>%
  group_by(age_group, risk_group) %>%
  mutate(prop_up_to_date_sim_period = risk_coverage_per_week/sum(risk_coverage_per_week),
         week = as.numeric(interval(as.Date('2024-06-29'), date) %/% weeks(1))) %>%
  filter(week != 0) %>%
  dplyr::select(week, date, age_group, risk_group, risk_coverage, prop_up_to_date_sim_period)


#write.csv(combined_with_risk %>% dplyr::select(-date), "data/clean-data/vaccine-coverage-by-age-and-risk-forModelValidation.csv")

plot_data <- melt(combined_with_risk, id = c("age_group", 'risk_group',  "date"))

plot_data %>% filter(variable == "risk_coverage", age_group == "18-29 years") %>%
  ggplot(aes(date, value, color = interaction(age_group, risk_group))) + 
  geom_line() + geom_point() +
  ylim(0,100) +
  xlab("Time") +
  ylab("% Uptake") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2024-06-28"), as.Date("2025-01-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("1st dose updated booster uptake per week over simulation period\nAge Group: 18-29 years")

