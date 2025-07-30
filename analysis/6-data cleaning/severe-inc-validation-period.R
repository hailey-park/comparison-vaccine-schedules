###################################################################################################
#Title: COVID-NET severe COVID-19 incidence cleaning 
#Author: Hailey Park
#Date: December 4, 2024
###################################################################################################

rm(list=ls())

#setwd(here::here())

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
library(data.table)

#Load in data
covidnet_char_data <- read.csv("raw-data/Patient_Characteristics_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20250321.csv")
covidnet_data <- read.csv("raw-data/Weekly_Rates_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20250322.csv")
prevalence_data <- read.csv("data/processed-data/risk_group_prevalence_estimates_1+comorbidity.csv")

#Estimating monthly incidence of COVID-19 hospitalizations by age group and risk, adjusting for COVID-19-reason hospitalizations
#NOTE:    Used data from CDC COVID-NET Surveillance system on COVID-19 Hospitalizations.
#         Link: https://data.cdc.gov/Public-Health-Surveillance/Weekly-Rates-of-Laboratory-Confirmed-COVID-19-Hosp/6jg4-xsqq/about_data
#         Link: https://data.cdc.gov/Public-Health-Surveillance/Patient-Characteristics-of-Laboratory-Confirmed-CO/bigw-pgk2/about_data 

#Pull relevant patient characteristics (reason admitted, underlying medical conditions)
covidnet_char_clean <- covidnet_char_data %>% 
  separate(Time, c("month", "year"), sep = "_") %>%
  filter(Time.Period == "Month",
         Sex == "Overall",
         Race.Ethnicity == "Overall",
         year %in% c("2023", "2024", "2025"),
         Estimate.Type == "Percent") %>%
  dplyr::select(Strata, Age.Category, COVID, month, year, Estimate) %>%
  mutate(month = case_when(month == "Novermber" ~ "November", 
                           month == "Febuary" ~ "February",
                           TRUE ~ month))

hosp_reason_data <- covidnet_char_clean %>% filter(Strata == "Admitted for COVID") %>% 
  rename(hosp_reason_estimate = Estimate,
         age_group = Age.Category) %>%
  dplyr::select(age_group, month, year, hosp_reason_estimate)

underlying_cond_data <- covidnet_char_clean %>% filter(Strata == "Any underlying condition")%>% 
  rename(underlying_cond_estimate = Estimate,
         age_group = Age.Category) %>%
  dplyr::select(age_group, month, year, underlying_cond_estimate)
    
combined_char_data <- merge(hosp_reason_data, underlying_cond_data, by = c("age_group", "month", "year"), all.x = TRUE)

#Split the 65+ age group into a 65-74 and 75+ age group
#below65 <- combined_char_data %>% filter(age_group %in% "≥65 years")
above65_74 <- combined_char_data %>% filter(age_group == "≥65 years") %>%
  mutate(age_group = "65-74 years")
above75 <- combined_char_data %>% filter(age_group == "≥65 years") %>%
  mutate(age_group = "≥75 years")

#Split the 18-49 year age group into 18-29 and 30-49
not_altered_age_groups <- combined_char_data %>% filter(!age_group %in% c("18-49 years", "≥65 years"))
age_18_29 <- combined_char_data %>% filter(age_group == "18-49 years") %>%
  mutate(age_group = "18-29 years")
age_30_49 <- combined_char_data %>% filter(age_group == "18-49 years") %>%
  mutate(age_group = "30-49 years")

combined_char_data <- as.data.frame(rbind(not_altered_age_groups, age_18_29, age_30_49, above65_74, above75))

#For missing data on % underlying medical conditions and % hospitalizations for COVID between November 2024 and January 2025, take average
average_data_2024 <- combined_char_data %>% filter(year == 2024, age_group != "Overall") %>%
  group_by(age_group) %>% summarise(hosp_reason_estimate = mean(hosp_reason_estimate),
                                    underlying_cond_estimate = mean(underlying_cond_estimate))

missing_data_months <- rbind(average_data_2024 %>% mutate(month = "November", year = 2024),
                             average_data_2024 %>% mutate(month = "December", year = 2024),
                             average_data_2024 %>% mutate(month = "January", year = 2025),
                             average_data_2024 %>% mutate(month = "February", year = 2025),
                             average_data_2024 %>% mutate(month = "March", year = 2025)) %>%
  dplyr::select(age_group, month, year, hosp_reason_estimate, underlying_cond_estimate)
  
combined_char_data <- as.data.frame(rbind(combined_char_data, missing_data_months))
                                    
#Clean COVID-NET severe inc data
severe_inc_clean <- covidnet_data %>% 
  mutate(week = as.Date(str_sub(X_WeekendDate, 1, 10)),
         weeks_since = as.numeric(interval(as.Date('2023-07-01'), week) %/% weeks(1)),
         month = format(week, '%B'),
         year = format(week, "%Y")) %>%
  filter(State == "COVID-NET",
         Season %in% c("2019-20","2020-21","2021-22", "2022-23", "2023-24", "2024-25"),
         AgeCategory_Legend %in% c("0-17 years (Children)", "18-49 years", "50-64 years", "65-74 years", "≥75 years"),
         Sex_Label == "All",
         Race_Label == "All",
         week >= as.Date("2023-07-01")) %>%
  dplyr::select(week, weeks_since, month, year, AgeCategory_Legend, WeeklyRate) %>%
  rename(age_group = AgeCategory_Legend,
         inc = WeeklyRate) 

#add extended age groups
age_18_29_inc <- severe_inc_clean %>% filter(age_group == "18-49 years") %>%
  mutate(age_group = "18-29 years")
age_30_49_inc <- severe_inc_clean %>% filter(age_group == "18-49 years") %>%
  mutate(age_group = "30-49 years")

severe_inc_clean <- data.frame(rbind(severe_inc_clean %>% filter(age_group != "18-49 years"), age_18_29_inc, age_30_49_inc))


#Clean prevalence data (3 risk groups)
prevalence_clean <- prevalence_data %>%
  pivot_wider(names_from = risk_group,
              values_from = prevalence) %>%
  mutate(age_group_prev = healthy + immunocompromised + `higher risk`)

#Add immunocompromised incidence data (https://www.cdc.gov/acip/downloads/slides-2024-10-23-24/03-COVID-Taylor-508.pdf)
immunocompromised_data <- data.frame(age_group = prevalence_clean$age_group,
                                     perc_immuno = c(14, 16, 16, 22, 15, 15))

#Combine severe inc data with patient characteristic info
combined_inc_with_char <- merge(merge(merge(severe_inc_clean, combined_char_data, by = c("age_group", "month", "year"), all.x = TRUE),
                                prevalence_clean, id = "age_group", all.x = TRUE),
                                immunocompromised_data, id = "age_group", all.x = TRUE)

adj_severe_inc_clean <- combined_inc_with_char %>%
  mutate(adj_inc = hosp_reason_estimate/100 * inc,
         higher_risk_inc = (underlying_cond_estimate - perc_immuno)/100 * adj_inc * age_group_prev/(`higher risk`),
         healthy_inc = (1 - underlying_cond_estimate/100) * adj_inc * age_group_prev/healthy,
         immunocompromised_inc = perc_immuno/100 * adj_inc * age_group_prev/(immunocompromised),
         age_group = case_when(age_group == "0-17 years (Children)" ~ "0-17 years",
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) %>%
  dplyr::select(age_group, week, weeks_since, inc, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc)


#write.csv(adj_severe_inc_clean %>% dplyr::select(week, weeks_since, age_group, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc), "data/clean-data/weekly-incidence-estimates-US-validationPeriod.csv")

plot_data <- melt(adj_severe_inc_clean, id = c("age_group", "week", "weeks_since"))


plot_data %>% 
  filter(variable %in% c("healthy_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  ylim(0,60) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2024-07-01"), as.Date("2025-01-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US\nData: COVID-NET\nRisk Group: Healthy")

###################################################################################################
actual_prevalence = entire_pop %>% group_by(age_group, risk_group) %>% summarise(prevalence = n()/10000000)

observed_severe_inc <- adj_severe_inc_clean %>% mutate(age_group = case_when(age_group == "≥75 years" ~ "75+ years",
                                                                            age_group == "0-17 years (Children)" ~ "0-17 years",
                                                                            TRUE ~ age_group))
  
overall_pop_inc <- merge(observed_severe_inc, actual_prevalence, by = "age_group", all.x = TRUE) %>%
  mutate(risk_adj_inc = case_when(risk_group == "healthy" ~ healthy_inc * prevalence,
                                  risk_group == "immunocompromised" ~ immunocompromised_inc * prevalence,
                                  risk_group == "higher risk" ~ higher_risk_inc * prevalence,
                                  TRUE ~ NA)) %>%
  group_by(week, weeks_since) %>% summarise(overall_inc = sum(risk_adj_inc), type = "Adjusted for Simulation Pop Distribution")

#estimate overall incidence from COVID-NET data
overall_covidnet_char_clean <- covidnet_char_data %>% 
  separate(Time, c("month", "year"), sep = "_") %>%
  filter(Time.Period == "Month",
         Sex == "Overall",
         Race.Ethnicity == "Overall",
         Age.Category == "Overall",
         year %in% c("2023", "2024"),
         Estimate.Type == "Percent") %>%
  dplyr::select(Strata, COVID, month, year, Estimate)

hosp_reason_data <- overall_covidnet_char_clean %>% filter(Strata == "Admitted for COVID") %>% 
  rename(hosp_reason_estimate = Estimate) %>%
  dplyr::select(month, year, hosp_reason_estimate)

underlying_cond_data <- overall_covidnet_char_clean %>% filter(Strata == "Any underlying condition")%>% 
  rename(underlying_cond_estimate = Estimate) %>%
  dplyr::select(month, year, underlying_cond_estimate)

combined_char_data <- merge(hosp_reason_data, underlying_cond_data, by = c("month", "year"), all.x = TRUE)

overall_severe_inc_clean <- covidnet_data %>% 
  mutate(week = as.Date(str_sub(X_WeekendDate, 1, 10)),
         weeks_since = as.numeric(interval(as.Date('2023-07-01'), week) %/% weeks(1)),
         month = format(week, '%B'),
         year = format(week, "%Y")) %>%
  filter(State == "COVID-NET",
         Season %in% c("2019-20","2020-21","2021-22", "2022-23", "2023-24"),
         AgeCategory_Legend %in% c("All"),
         Sex_Label == "All",
         Race_Label == "All",
         week >= as.Date("2023-07-01")) %>%
  dplyr::select(week, weeks_since, month, year, WeeklyRate) %>%
  rename(inc = WeeklyRate)

#Combine severe inc data with patient characteristic info
combined_inc_with_char <- merge(overall_severe_inc_clean, combined_char_data, by = c("month", "year"), all.x = TRUE)
 
adj_severe_inc_clean <- combined_inc_with_char %>%
  mutate(overall_inc = hosp_reason_estimate/100 * inc,
         type = "COVID-NET Overall Pop Estimate") %>%
  dplyr::select(week, weeks_since, overall_inc, type)

combined <- rbind(overall_pop_inc, adj_severe_inc_clean)

ggplot(data = combined, aes(x = as.Date(week), y = overall_inc, color = type)) + geom_line() +
  ylim(0, 8) + ylab("Weekly Severe Incidence (per 100,000)") +
  xlab("Time") + xlim(as.Date("2023-07-01"), as.Date("2024-07-01")) +
  ggtitle("How Simulated Population Distribution Affects Overall Severe Incidence Estimate")



covid_net_overall <- read.csv("data/clean-data/weekly-incidence-estimates-US-overall.csv")[,-1]
