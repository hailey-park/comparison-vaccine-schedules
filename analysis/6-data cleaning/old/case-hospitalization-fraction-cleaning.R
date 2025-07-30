###################################################################################################
#Title: Case-hospitalization fraction cleaning 
#Author: Hailey Park
#Date: December 18, 2024
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
cases_deaths_cdc <- read.csv("data/raw-data/Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Updated__Bivalent__Booster_Status_20241218.csv")


#Clean data
case_death_frac <- melt(cases_deaths_cdc %>%
  filter(month == "MAR 2023", age_group != "all_ages") %>%
  mutate(revised_age_group = case_when(age_group %in% c("0.5-4", "5-11", "12-17") ~ "0-17 years",
                                       age_group == "18-29" ~ "18-29 years",
                                       age_group == "30-49" ~ "30-49 years",
                                       age_group == "50-64" ~ "50-64 years",
                                       age_group == "65-79" ~ "65-74 years",
                                       age_group == "80+" ~ "75+ years",
                                       TRUE ~ NA),
         count = vaccinated_with_outcome + unvaccinated_with_outcome,
         pop = vaccinated_population + unvaccinated_population) %>%
  group_by(outcome, mmwr_week, revised_age_group) %>%
  summarise(across(c(count, pop), list(sum))) %>%
  group_by(outcome, revised_age_group) %>%
  summarise(prop_outcome = mean(count_1/pop_1 * 100000)) %>%
  pivot_wider(names_from = revised_age_group,
              values_from = prop_outcome) %>%
  mutate(`75+ years` = 0.15 * `65-74 years` + 0.85 * `75+ years`), id = "outcome") %>%
  rename(revised_age_group = variable, prop_outcome = value) %>%
  pivot_wider(names_from = outcome,
              values_from = prop_outcome) %>%
  add_column(hosp = c(0.8, 1.2, 1.6, 4.45, 10.07, 31.6)) %>%
  mutate(case_death_frac = (death/case),
         case_hosp_frac = (hosp/case),
         cases_per_1_death = 1/case_death_frac,
         cases_per_1_hosp = 1/case_hosp_frac)

