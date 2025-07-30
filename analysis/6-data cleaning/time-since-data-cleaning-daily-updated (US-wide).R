###################################################################################################
#Title: Cleaning Case Data and Vaccine Data Used for 'Time Since Last' Estimation
#Author: Hailey Park
#Date: December 12, 2024
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
#Case data from this website (https://covidestim.org/))
#Vaccine data through May 10, 2023 from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/about_data)
#Vaccine data from June 14 - October 11 2023 from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccines-Up-to-Date-Status/9b5z-wnve/about_data)

cases <- read.csv("data/raw-data/state - from covidestim/estimates.csv")
vaccines_may2023 <- read.csv("data/raw-data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional_20241217.csv")
vaccines_oct2023 <- read.csv("data/raw-data/COVID-19_Vaccines_Up_to_Date_Status_20240702.csv")

###################################################################################################
#Cases by Day
cases_by_day <- cases %>% mutate(date = as.Date(date)) %>% 
  filter(date >= as.Date('2022-01-01') & date <= as.Date('2023-07-01')) %>%
  group_by(date) %>% summarise(total_cases = sum(fitted_cases)) %>% 
  complete(date = full_seq(date, 1)) %>% 
  fill(c(total_cases)) %>% 
  mutate(perc_cases = total_cases/sum(total_cases))


ggplot(data = cases_by_day, aes(x = date, y = total_cases)) +
  geom_point() + geom_line()
################################################################################################### 
#Vaccine Doses by Day (for primary series and 1st booster coverage, 2nd booster coverage will be calculated later)
doses_by_day <- vaccines_may2023 %>% 
  filter(Location %in% c("US"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y%H")) %>%
  filter(date >= as.Date('2020-12-15') & date <= as.Date('2023-07-01')) %>%
  group_by(age_group, date) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                          total_primary_series = sum(Series_Complete_Yes, na.rm=TRUE),
                                          total_boosted_1st = sum(Booster_Doses, na.rm=TRUE)) %>%
  group_by(age_group) %>% mutate(primary_series_today = c(total_primary_series[1], (total_primary_series - lag(total_primary_series, 1))[-1]),
                                 boosted_today = c(total_boosted_1st[1], (total_boosted_1st - lag(total_boosted_1st, 1))[-1]))

#Estimating % of people who completed their 1st booster dose at each day, from December 15, 2020 - May 10, 2023
#We are assuming no 1st booster vaccination happened between May 10-July 1, 2023 since data is unavailable.
booster_doses_by_day <- doses_by_day %>% filter(date >= as.Date("2021-09-22")) %>%
  mutate(day = as.Date(date), week = floor_date(as.Date(date), unit = "week")) %>% group_by(week, age_group) %>% mutate(boosted_today = mean(boosted_today)) %>%
  group_by(day, age_group) %>% summarise(total_vaccinated = sum(boosted_today)) %>%
  group_by(age_group) %>% mutate(perc_doses = total_vaccinated/sum(total_vaccinated)) %>%
  dplyr::select(day, age_group, perc_doses) %>%
  mutate(age_group = if_else(age_group == "18-49 years", "18-29 years", age_group))
  
age_30_49_only <- booster_doses_by_day %>% filter(age_group == "18-29 years") %>% mutate(age_group = "30-49 years")
booster_doses_by_day_final <- data.frame(rbind(booster_doses_by_day, age_30_49_only))

# Estimating the % of people who completed primary series at each day, from December 15, 2020 - May 10, 2023. 
# We are assuming no primary series vaccination happened between May 10, 2023 - July 1, 2023. In reality, very
# little actually happened so fine assumption to make.
fully_vax_doses_by_day <- doses_by_day %>% 
  mutate(day = as.Date(date), week = floor_date(as.Date(date), unit = "week")) %>% group_by(week, age_group) %>% mutate(primary_series_today = mean(primary_series_today)) %>%
  group_by(day, age_group) %>% summarise(total_vaccinated = sum(primary_series_today)) %>%
  group_by(age_group) %>% mutate(perc_doses = total_vaccinated/sum(total_vaccinated)) %>%
  dplyr::select(day, age_group, perc_doses) %>%
  mutate(age_group = if_else(age_group == "18-49 years", "18-29 years", age_group))

age_30_49_only <- fully_vax_doses_by_day %>% filter(age_group == "18-29 years") %>% mutate(age_group = "30-49 years")
fully_vax_doses_by_day_final <- data.frame(rbind(fully_vax_doses_by_day, age_30_49_only))

################################################################################################### 
#Vaccine doses by day through May 2023 (2nd boosters monovalent)
doses_by_day_may2023 <- vaccines_may2023 %>% 
  filter(Location %in% c("US"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y%H"),
         day = as.Date(date),
         week = floor_date(as.Date(date), unit = "week")) %>%
  filter(date >= as.Date("2021-12-25") & date <= as.Date('2024-06-01')) %>%
  group_by(age_group, date, day, week) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                                total_boosted_1st = sum(Booster_Doses, na.rm=TRUE),
                                                total_boosted_2nd = sum(Second_Booster, na.rm=TRUE)) %>%
  # group_by(age_group, week) %>% mutate(total_boosted_1st = mean(total_boosted_1st),
  #                                      total_boosted_2nd = mean(total_boosted_2nd)) %>%
  group_by(age_group, day) %>% summarise(total_pop = mean(total_pop, na.rm=TRUE),
                                        total_boosted_1st = max(total_boosted_1st),
                                        total_boosted_2nd = max(total_boosted_2nd)) %>%
  group_by(age_group) %>% mutate(day_num = as.numeric(difftime(day, as.Date("2021-12-26"), units = "days"))) %>%
  filter(day != as.Date("2021-12-19"))

#Vaccine doses by month from June 2023 through October 2023 (up-to-date boosters)
doses_by_month_oct2023 <- vaccines_oct2023 %>% 
  filter(!Location %in% c("PW", "PR", "GU", "VI", "TT", "MP", "AS", "MH"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(date >= as.Date('2022-01-01') & date <= as.Date('2023-07-12')) %>%
  group_by(age_group, date) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                          total_up_to_date = sum(Up_to_date, na.rm=TRUE))

#Imputing the missing data on 2nd booster uptake for the 0-17 years and 18-49 years age group

  #Using the 50-64 year age group data to impute the 2nd booster uptake trend over time for younger age groups
  #Using 1st booster relationship @ July 2023 to fix relationship for 2nd booster
age_ratio_1st_booster <- doses_by_day_may2023 %>% filter(day == as.Date("2023-05-10"), !age_group %in% c("65-74 years", "75+ years")) %>%
  mutate(ratio = (total_boosted_1st/total_pop)/(31044427/63659835)) %>%
  dplyr::select(age_group, ratio)

imputed_younger_age_coverage_data <- melt(doses_by_day_may2023 %>% mutate(coverage = (total_boosted_2nd)/(total_pop)) %>%
  filter(!age_group %in% c("65-74 years", "75+ years")) %>%
  dplyr::select(age_group, day, coverage) %>% 
  pivot_wider(names_from = age_group,
              values_from = coverage) %>%
  rowwise() %>%
  mutate(`0-17 years` = max(`50-64 years` * (age_ratio_1st_booster %>% filter(age_group == "0-17 years"))$ratio - 0.0080284, 0),
         `18-49 years` = max(`50-64 years` * (age_ratio_1st_booster %>% filter(age_group == "18-49 years"))$ratio - 0.04608035, 0)),
  id = "day") %>% rename(age_group = variable, coverage = value)

#Imputing the missing data on 2nd booster uptake between June-July of 2023 (all age groups)
#Taking the May 2023 estimate from the 2nd boosters CDC dataset and the July 2023 estimate from
#the up-to-date coverage estimate, and assuming linear increase in coverage between those dates.
missing_data <- data.frame(day = rep(seq(as.Date("2023-05-11"), as.Date("2023-07-01"), by = "days"), 5),
                           age_group = rep(c("0-17 years", "18-49 years", "50-64 years", "65-74 years", "75+ years"), each = 52),
                           coverage = c(c(0.037 + c(1:52) * (0.0451 - 0.037)/52), #0-17 years
                                        c(0.0845 + c(1:52) * (0.101 - 0.0845)/52), #18-49 years
                                        rep(0.2144, 52), #50-64 years
                                        rep(0.411, 52), #65-74 years
                                        rep(0.436, 52))) #75+ years

#Combine imputed data together
combined_coverage_data <- merge(merge(doses_by_day_may2023, imputed_younger_age_coverage_data, by = c("age_group", "day"), all.x = TRUE),
                                missing_data, by = c("age_group", "day"), all = TRUE) %>%
  rowwise() %>% mutate(coverage = max(c(coverage.x, coverage.y), na.rm = TRUE),
                       coverage = if_else(coverage == -Inf, total_boosted_2nd/total_pop, coverage)) 
  
#Estimating the % of people who get a 2nd booster each day
booster_doses_2nd_by_day <- combined_coverage_data %>%
  group_by(age_group) %>%
  mutate(coverage_today = c(coverage[1], (coverage - lag(coverage, 1))[-1])) %>%
  rowwise() %>% mutate(coverage_today = max(coverage_today, 0)) %>%
  mutate(week = floor_date(as.Date(day), unit = "week")) %>% group_by(week, age_group) %>% mutate(coverage_today = mean(coverage_today)) %>%
  ungroup() %>% group_by(age_group) %>% mutate(perc_doses = coverage_today/sum(coverage_today)) %>% 
  dplyr::select(day, age_group, perc_doses) %>%
  mutate(age_group = if_else(age_group == "18-49 years", "18-29 years", age_group)) %>%
  filter(day >= as.Date("2022-01-01"))

age_30_49_only <- booster_doses_2nd_by_day %>% filter(age_group == "18-29 years") %>% mutate(age_group = "30-49 years")
booster_doses_2nd_by_day_final <- data.frame(rbind(booster_doses_2nd_by_day, age_30_49_only))


################################################################################################### 
#write as .csv
# write.csv(cases_by_day, "data/clean-data/cases_by_day.csv")
# write.csv(booster_doses_by_day_final, "data/clean-data/booster_1st_doses_by_day.csv")
# write.csv(booster_doses_2nd_by_day_final, "data/clean-data/booster_2nd_doses_by_day_updated.csv")
# write.csv(fully_vax_doses_by_day_final, "data/clean-data/fully_vax_doses_by_day.csv")

#Plot coverages over time
ggplot(data = booster_doses_2nd_by_day_final %>% filter(age_group == "0-17 years"), aes(x = day, y = perc_doses, color = age_group)) +
  geom_line() +
  ylim(0, NA) +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2021-09-01"), as.Date("2023-07-01")),
               breaks = "3 months") +
  xlab("Time") +
  ylab("Proportion") +
  ggtitle("Proportion of Completed 1st Booster Over Time")

ggplot(data = combined_coverage_data, aes(x = day, y = coverage, color = age_group)) +
  geom_line() +
  ylim(0, 1) +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2021-09-01"), as.Date("2023-07-01")),
               breaks = "3 months") +
  xlab("Time") +
  ylab("Proportion") +
  ggtitle("Proportion of Completed 2nd Booster Over Time")
  



