###################################################################################################
#Title: Cleaning Vaccine Coverage Data Used
#Author: Hailey Park
#Date: June 19, 2024
###################################################################################################

rm(list=ls())

#setwd(here::here())

#Load libraries
library(tidyverse)
library(reshape2)

#Read in data
  #Vaccine data from June 14 - October 11 2023 from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccines-Up-to-Date-Status/9b5z-wnve/about_data)
  #Vaccine data from  October 2023 - December 2024 from this website (https://data.cdc.gov/Vaccinations/Weekly-Cumulative-COVID-19-Vaccination-Coverage-an)
  #Vaccine data for children from  October 2023 - December 2024 from this website (https://data.cdc.gov/Child-Vaccinations/Weekly-Parental-Intent-for-Vaccination-and-Cumulat/ker6-gs6z/about_data)
vaccines_oct2023 <- read.csv("raw-data/COVID-19_Vaccines_Up_to_Date_Status_20240702.csv")
vaccine_dec2024 <- read.csv("raw-data/Weekly_Cumulative_COVID-19_Vaccination_Coverage_and_Intent__Overall__by_Selected_Demographics_and_Jurisdiction__Among_Adults_18_Years_and_Older_20250107.csv")
vaccine_dec2024_child <- read.csv("raw-data/Weekly_Parental_Intent_for_Vaccination_and_Cumulative_Percentage_of_Children_6_Months_-17_Years_Who_are_Up_to_date_with_the_COVID-19_Vaccines_by_Season__United_States_20250107.csv")

  #Vaccine data from October 2023 - July 2024 from this website (https://www.cdc.gov/covidvaxview/interactive/adults.html?CDC_AAref_Val=https://www.cdc.gov/vaccines/imz-managers/coverage/covidvaxview/interactive/adults.html)
vaccines_risk <- read.csv("raw-data/National_Immunization_Survey_Adult_COVID_Module__NIS-ACM___COVIDVaxViews__Data___Centers_for_Disease_Control_and_Prevention__cdc.gov__20250114.csv")

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
vaccines_by_risk <- vaccines_risk %>% filter(Geography.Type == "National Estimates", 
                                             Group.Name == "Health condition associated with higher risk for COVID-19 (any)",
                                             Indicator.Name == "Vaccination Coverage and Intent",
                                             Indicator.Category == "Received a 2023-2024 COVID-19 vaccine dose (among all adults 18+)") %>%
  mutate(month = case_when(Time.Period == "April 1 - April 27" ~ as.Date("2024-04-01"),
                           Time.Period == "April 28 - May 25" ~ as.Date('2024-05-01'),
                           Time.Period == "December 31 - January 27" ~ as.Date('2024-01-01'),
                           Time.Period == "February 25 - March 31" ~ as.Date('2024-03-01'),
                           Time.Period == "January 28 - February 24" ~ as.Date("2024-02-01"),
                           Time.Period == "July 1 - July 27" ~ as.Date('2024-07-01'),
                           Time.Period == "May 26 - June 30" ~ as.Date('2024-06-01'),
                           Time.Period == "November 26 - December 30" ~ as.Date('2023-12-01'),
                           Time.Period == "October 29 - November 25" ~ as.Date('2023-11-01'),
                           Time.Period == "September 24 - October 28" ~ as.Date("2023-10-01"),
                           TRUE ~ NA)) %>%
  dplyr::select(c(Group.Category, month, Estimate....)) %>% rename(high_risk = Group.Category, vaccine_coverage = Estimate....) %>%
  pivot_wider(names_from = high_risk,
              values_from = vaccine_coverage) %>%
  mutate(ratio = Yes/No)

average_ratio <- mean(vaccines_by_risk$ratio)

#Vaccine doses by month from June 2023 through October 2023
max_min_perc <- vaccines_oct2023 %>% 
  filter(Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs") ~ "18-29 years",
                               Demographic_Category %in% c("Ages_25-49_yrs") ~ "30-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(date >= as.Date('2023-07-01') & date <= as.Date('2023-10-01')) %>%
  group_by(age_group, date) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                                total_up_to_date = sum(Up_to_date, na.rm=TRUE)) %>%
  group_by(age_group) %>% summarise(max_perc = max(total_up_to_date/total_pop),
                                    min_perc = min(total_up_to_date/total_pop))

doses_by_month_oct2023 <- data.frame(date = rep(seq(as.Date("2023-07-01"), as.Date("2023-09-29"), by = 'week'), times = 6),
                                     age_group = rep(c("0-17 years", "18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"), each = 13),
                                     perc_up_to_date_week = rep(((max_min_perc$max_perc - max_min_perc$min_perc)/(13)) * 100, each = 13)) %>%
  group_by(age_group) %>% mutate(coverage = cumsum(perc_up_to_date_week))

#Vaccine doses by month from October 2023 through December 2024 (ADULTS)
doses_by_month_dec2024_adult <- vaccine_dec2024 %>% filter(Geographic.Level == "National", Demographic.Level == "Age", indicator_label == "Up-to-date", Demographic.Name %in% c("18-29 years", "30-39 years", "40-49 years", "50-64 years", "65-74 years", "75+ years")) %>%
  mutate(date = as.Date(Week_ending, format = "%m/%d/%Y"),
         age_group = if_else(Demographic.Name %in% c("30-39 years", "40-49 years"), "30-49 years", Demographic.Name)) %>% group_by(date, age_group) %>%
  summarise(perc = mean(Estimate)) %>%
  group_by(age_group) %>% mutate(perc_up_to_date_week = c((perc - lag(perc, 1))[2], (perc - lag(perc, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, perc, perc_up_to_date_week)) %>%
  filter(date < as.Date("2024-07-01")) %>% dplyr::select(date, age_group, perc_up_to_date_week) %>%
  group_by(age_group) %>% mutate(coverage = cumsum(perc_up_to_date_week))
  
#Vaccine doses by month from October 2023 through December 2024 (CHILDREN)
doses_by_month_dec2024_child <- vaccine_dec2024_child %>% filter(Geographic.Level == "National", Demographic_Level == "Age", Indicator_label == "Up-to-date", Demographic.Name %in% c("6 months-4 years", "5-17 years")) %>%
  mutate(date = as.Date(Week_ending)) %>% group_by(date) %>%
  summarise(perc = weighted.mean(Estimate, Unweighted.Sample.Size) * 100) %>%
  mutate(perc_up_to_date_week = c((perc - lag(perc, 1))[2], (perc - lag(perc, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, perc, perc_up_to_date_week),
         age_group = "0-17 years") %>%
  filter(date < as.Date("2024-07-01")) %>% dplyr::select(date, age_group, perc_up_to_date_week) %>%
  group_by(age_group) %>% mutate(coverage = cumsum(perc_up_to_date_week))

#Combine vaccine uptake data
combined <- rbind(doses_by_month_oct2023, doses_by_month_dec2024_adult, doses_by_month_dec2024_child)

#ALTERNATIVE UPTAKE SCENARIOS (pessimistic)
# 0-17 years: 0%
# 18-49 years: 5%
# 50-64 years: 10%
# 65+ years: 20%
combined_pessimistic <- combined %>%
  mutate(coverage = case_when(age_group == "0-17 years" ~ 0,
                              age_group == "18-29 years" ~ coverage * 5/11.3,
                              age_group == "30-49 years" ~ coverage * 5/15.95,
                              age_group == "50-64 years" ~ coverage * 10/21.7,
                              age_group == "65-74 years" ~ coverage * 20/37.9,
                              age_group == "75+ years" ~ coverage * 20/37.6)) %>%
  group_by(age_group) %>% mutate(perc_up_to_date_week = c((coverage - lag(coverage, 1))[2], (coverage - lag(coverage, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, coverage, perc_up_to_date_week))

#ALTERNATIVE UPTAKE SCENARIOS (optimistic)
# 0-17 years: 50%
# 18-49 years: 60%
# 50-64 years: 70%
# 65+ years: 80%
combined_optimistic <- combined %>%
  mutate(coverage = case_when(age_group == "0-17 years" ~ coverage * 50/14.3,
                              age_group == "18-29 years" ~ coverage * 60/11.3,
                              age_group == "30-49 years" ~ coverage * 60/15.95,
                              age_group == "50-64 years" ~ coverage * 70/21.7,
                              age_group == "65-74 years" ~ coverage * 80/37.9,
                              age_group == "75+ years" ~ coverage * 80/37.6)) %>%
  group_by(age_group) %>% mutate(perc_up_to_date_week = c((coverage - lag(coverage, 1))[2], (coverage - lag(coverage, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, coverage, perc_up_to_date_week))

#Calculate healthy group vaccine uptake (CHANGE `combined` with `combined_pessimistic` or `combined_optimistic`)
healthy_vax_coverage <- merge(combined %>% mutate(risk_group = "healthy"),
                              prevalence_adj, by = "age_group", all.x = TRUE) %>%
  mutate(risk_coverage = coverage/(Yes * 1.8 + No),
         risk_coverage_per_week = perc_up_to_date_week/(Yes * 1.8 + No))

#Add risk group information
combined_with_risk <- rbind(healthy_vax_coverage,
                            healthy_vax_coverage %>% mutate(risk_group = "immunocompromised",
                                                            risk_coverage = risk_coverage * 1.8,
                                                            risk_coverage_per_week = risk_coverage_per_week * 1.8),
                            healthy_vax_coverage %>% mutate(risk_group = "higher risk",
                                                            risk_coverage = risk_coverage * 1.8,
                                                            risk_coverage_per_week = risk_coverage_per_week * 1.8)) %>%
  mutate(week = as.numeric(interval(as.Date('2023-07-01'), date) %/% weeks(1))) %>%
  filter(week != 0) %>%
  group_by(age_group, risk_group) %>%
  mutate(prop_up_to_date_sim_period = risk_coverage_per_week/sum(risk_coverage_per_week)) %>%
  dplyr::select(week, date, age_group, risk_group, risk_coverage, prop_up_to_date_sim_period)


#write.csv(combined_with_risk, "data/clean-data/vax-uptake-scenarios/optimistic-vaccine-coverage-by-age-and-risk.csv")

plot_data <- melt(combined_with_risk, id = c("age_group", 'risk_group',  "date"))

plot_data %>% filter(variable == "prop_up_to_date_sim_period", age_group == "18-29 years") %>%
  ggplot(aes(date, value, color = interaction(age_group, risk_group))) + 
  geom_line() + geom_point() +
  ylim(0,NA) +
  xlab("Time") +
  ylab("% Uptake") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2024-07-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("1st dose updated booster uptake per week over simulation period\nAge Group: 18-29 years")

###################################################################################################
#2nd dose uptake (see screenshot in timed-booster slack from nathan)
second_dose_annually <- data.frame(date = rep(seq(as.Date('2024-03-02'), as.Date("2024-07-01"), by = "weeks"), 6 * 3),
                                   age_group = rep(c("0-17 years", "18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"), each = 18, times = 3),
                                   risk_group = rep(c("healthy", "immunocompromised", "higher risk"), each = 6 * 18)) %>%
  group_by(age_group, risk_group) %>% mutate(coverage = case_when(age_group %in% c("65-74 years", "75+ years") ~ c(cumsum(rep(3.8/4,4)), 3.8 + cumsum(rep((5.7 - 3.8)/5,5)), 5.7+cumsum(rep((7.3 - 5.7)/4,4)), 7.3+cumsum(rep((8.9 - 7.3)/5,5))),
                                                                  age_group == "0-17 years" ~ 0,
                                                                  risk_group != "healthy" ~ c(cumsum(rep(2.5/4,4)), 2.5 + cumsum(rep((3.3 - 2.5)/5,5)), 3.3 + cumsum(rep((4.5 - 3.3)/4,4)), 4.5 + cumsum(rep((5.4 - 4.5)/5,5))),
                                                                  TRUE ~ 0),
                                             perc_up_to_date_week = c((coverage - lag(coverage, 1))[2], (coverage - lag(coverage, 1))[-1])) %>% rowwise() %>%
  mutate(perc_up_to_date_week = if_else(perc_up_to_date_week < 0, coverage, perc_up_to_date_week)) %>%
  group_by(age_group, risk_group) %>%
  #rowwise() %>%
  mutate(prop_up_to_date_sim_period = if_else(is.na(perc_up_to_date_week/sum(perc_up_to_date_week)), 0, perc_up_to_date_week/sum(perc_up_to_date_week)),
         week = as.numeric(interval(as.Date('2023-07-01'), date) %/% weeks(1))) %>%
  filter(week != 0) %>%
  dplyr::select(week, date, age_group, risk_group, coverage, prop_up_to_date_sim_period)


#ALTERNATIVE UPTAKE SCENARIOS (pessimistic)
# Higher risk/immunocompromised and <65 years: 2.7% (realistic is 5.4%)
# 65+ years: 4.5% (realistic is 8.9%)
combined_pessimistic <- second_dose_annually %>%
  mutate(coverage = case_when(age_group %in% c("65-74 years", "75+ years") ~ coverage * (4.5/8.9),
                              age_group == "0-17 years" ~ 0,
                              risk_group != "healthy" ~ coverage * 2.7/5.4,
                              TRUE ~ 0)) %>%
  group_by(age_group, risk_group) %>%
  mutate(perc_up_to_date_week = c((coverage - lag(coverage, 1))[2], (coverage - lag(coverage, 1))[-1])) 

#ALTERNATIVE UPTAKE SCENARIOS (optimistic)
# Higher risk/immunocompromised and <65 years: 60-70% (realistic is 5.4%)
# 65+ years: 80% (realistic is 8.9%)
combined_optimistic <- second_dose_annually %>%
  mutate(coverage = case_when(age_group %in% c("65-74 years", "75+ years") ~ coverage * (80/8.9),
                              age_group == "0-17 years" ~ 0,
                              age_group %in% c("18-29 years", "30-49 years") & risk_group != "healthy" ~ coverage * 60/5.4,
                              age_group %in% c("50-64 years") & risk_group != "healthy" ~ coverage * 70/5.4,
                              TRUE ~ 0)) %>%
  group_by(age_group, risk_group) %>%
  mutate(perc_up_to_date_week = c((coverage - lag(coverage, 1))[2], (coverage - lag(coverage, 1))[-1])) 



#write.csv(combined_optimistic, "data/clean-data/vax-uptake-scenarios/optimistic-2nd-dose-annually-vaccine-coverage-by-age-and-risk.csv")                                   

plot_data <- melt(combined_optimistic, id = c("age_group", 'risk_group',  "date"))

plot_data %>% filter(variable == "coverage", age_group == "0-17 years") %>%
  ggplot(aes(date, value, color = interaction(age_group, risk_group))) + 
  geom_line() + geom_point() +
  ylim(0,100) +
  xlab("Time") +
  ylab("% Uptake") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2024-07-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("2nd dose updated booster uptake per week over simulation period\nAge Group: 0-17 years")


first_dose_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/optimistic-vaccine-coverage-by-age-and-risk.csv")[,-1] %>% mutate(dose = "1st dose") %>% rename(coverage = risk_coverage)
second_dose_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/optimistic-2nd-dose-annually-vaccine-coverage-by-age-and-risk.csv")[,-1] %>% mutate(dose = "2nd dose") %>% dplyr::select(-perc_up_to_date_week)

combined <- rbind(first_dose_coverage, second_dose_coverage)

plot_data <- melt(combined, id = c("age_group", 'risk_group', 'dose', "date"))

plot_data %>% filter(variable == "coverage", age_group == "75+ years") %>%
  ggplot(aes(as.Date(date), value, color = interaction(age_group, risk_group), shape = dose)) + 
  geom_line() + geom_point() +
  ylim(0,100) +
  xlab("Time") +
  ylab("% Uptake") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2024-07-01")),
               breaks = "1 months") +
  labs(color = "Age Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Realistic booster uptake per week over simulation period\nAge Group: 0-17 years")


