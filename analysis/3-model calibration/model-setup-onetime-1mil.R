###################################################################################################
#Title: Model setup (load one-time)
#Author: Hailey Park
#Date: March 8, 2025
###################################################################################################

#Historical vaccine coverage by age and risk group. Calculated proportion of population who receive vaccine each week of simulation period 
vaccine_coverage_by_age <- read.csv("data/clean-data/vax-uptake-scenarios/realistic-vaccine-coverage-by-age-and-risk.csv")[,-1] 
second_dose_annually_vaccine_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/realistic-2nd-dose-annually-vaccine-coverage-by-age-and-risk.csv")[,-1]
vaccine_coverage_by_age_validation <- read.csv("data/clean-data/vaccine-coverage-by-age-and-risk-forModelValidation.csv")[,-1]

#Severe incidence estimates by age (COVID-NET data, cleaned)
weekly_severe_incidence <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod.csv")[,-1] %>%
  mutate(age_group = case_when(age_group == "0-17 years (Children)" ~  "0-17 years", 
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) 

#Ratios between age-specific case hospitalization fractions
# NOTE: 1. The baseline case-hospitalization fraction `case_hosp_frac` is a calibrated parameter that sets the relationship between
#          severe and non-severe cases.
#       2. The case-hospitalzrelationship between all the The reference 50-64 years
age_ratio_case_hosp_frac <- data.frame(age_group = c("0-17 years", "18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),
                                       case_hosp_frac = c(0.155, 0.4283, 1.545, 5.9, 13, 30)) %>%
  mutate(ratio = case_hosp_frac/5.9)

#Contact matrix (POLYMOD)
contact_matrix <- read.csv("data/processed-data/contact matrix updated.csv") 

#Waning curves
waning_data_clean <- setDT(read.csv("data/clean-data/waning_data_clean.csv")[,-1])
waning_data_clean_alt <- setDT(read.csv("data/clean-data/waning_data_clean_alt.csv")[,-1]) %>%
  filter(weeks %% 2 == 0) %>% mutate(weeks = weeks/2)
setkeyv(waning_data_clean_alt, c("immuno", "weeks", "prior_inf", "prior_vacc"))

#Calculating total current infections at model initialization, based on average severe and nonsevere incidence estimates at July 2023
average_severe_incidence <- melt(weekly_severe_incidence %>% filter(weeks_since < 4) %>% group_by(age_group) %>% 
  summarise(across(healthy_inc:higher_risk_inc, mean)), id = "age_group") %>% 
  rowwise() %>% dplyr::mutate(risk_group = strsplit(as.character(variable), split="_")[[1]][1],
                              risk_group = if_else(risk_group == "higher", "higher risk", risk_group)) %>%
  rename(severe_inc = value) %>% dplyr::select(-variable)

#Read in initialized population file
entire_pop <-  read.csv(paste0("data/clean-data/entire_population_model_initialization_daily_1mil_updated.csv"))[,-1] 

#Read vaccine-uptake-assignment script and run vaccine assignment based on historical coverage on entire population
source(here::here("vaccine-uptake-assignment-calibration.R"))

entire_pop_with_vax_assigned <- realistic_vax_assignment(entire_pop)

############################################################################################