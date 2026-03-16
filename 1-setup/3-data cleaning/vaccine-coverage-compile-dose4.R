########################################################################################################################
# Title: Make dose4_data file
# Author: Kate Bubar
# Date: March 13, 2026
# 
# Overview: Create file storing dose4_data, assuming same timing and coverage as dose2_coverage with 
#           updated dates and shorter time span
#
# Details:
#           Dose 2 spanned 4 months (March-June 2024), sim days 240-365
#           Dose 4 will span 2 months (March-April 2025), sim days 609-665
#           Change first week from 35 -> 87, 03-02-2024 -> 03-01-2025
#                  last week is week 95, 04-26-2025
#########################################################################################################################

library(dplyr)
library(lubridate)

dose2_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/realistic-2nd-dose-annually-vaccine-coverage-by-age-and-risk.csv")[,-1]

dose4_coverage <- dose2_coverage %>%
  mutate(
    week = week + 52,
    date = as.character(as.Date(date, format = "%Y-%m-%d") + 364)
  ) %>%
  filter(week <= 95)

#write.csv(dose4_coverage, "data/clean-data/vax-uptake-scenarios/realistic-dose4-coverage-by-age-and-risk.csv")
