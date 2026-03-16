###################################################################################################
# Title: Run Model - 1 iteration, written to run locally
# Author: Kate Bubar
# Date: January 29, 2026
###################################################################################################

rm(list=ls())

#start <- Sys.time()
#Load libraries
library(tidyverse)
library(data.table)

###################################################################################################
#LOAD in scripts (only needs to be done once)
setwd("~/comparison-vaccine-schedules/")
source("functions/vaccine-uptake-assignment-historical.R")
source("functions/simulation-functions-gridsearch.R") # includes how long to run model
#source(here::here("simulation-functions-updated-t0.R")) #FIXME: KMB add details of when to run this or what is output

###################################################################################################
#READ in parameters and model inputs
# FIXME: KMB - add code to choose pop size
start.time <- Sys.time()
#source("config/model_parameters.R")

# Burn-in period
lambda_1 <- 0.49
lambda_2 <- 0.42
lambda_3 <- 0.45
lambda_4 <- 0.46
lambda_5 <- 0.59
lambda_6 <- 0.66
lambda_7 <- 0.58
lambda_8 <- 0.62
lambda_9 <- 0.57

# Phase 1
lambda_10 <- 0.52#-0.05
lambda_11 <- 0.52#-0.05
lambda_12 <- 0.52#-0.05
lambda_13 <- 0.52#-0.05 # changed from 0.55!
lambda_14 <- 0.52#-0.05
lambda_15 <- 0.52#-0.05

# Phase 2
lambda_16 <- 0.55
lambda_17 <- 0.6
lambda_18 <- 0.65

lambda_19 <- 0.6
lambda_20 <- 0.5
lambda_21 <- 0.5
lambda_22 <- 0.5
baseline_case_hosp_frac <- 0.0115 # 0.012

source("functions/simulation-inputs-historical-5mil.R") # 2.5 min
print(paste0("time: ", round(Sys.time() - start.time, 2)))

#ASSIGN vaccine uptake and timing to each individual in clean_df
start.time <- Sys.time()
set.seed(1)
df_sim <- historical_vax_assignment(clean_df[[1]], dose1_coverage, dose2_coverage, dose3_coverage, dose4_coverage) # 1 min


#RUN the simulation
#   NOTE: vax distribution used in semiannual_strat_8_9_10 corresponds to the observed historical distributions
results <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = 274, rel_ve_check = FALSE)
print(paste0("time: ", round(Sys.time() - start.time, 2)))

1-(results$severe_vax1/results$total_vax1)/(results$severe_no_vax1/results$total_no_vax1)


###################################################################################################
#SUMMARIZE the results to output results aggregated by age group
#   NOTE: If you want results by age and risk, output results as is
age_group_results <- results %>%
  group_by(age_group) %>%
  summarise(across(-risk_group, sum)) %>%
  pivot_longer(-age_group, names_to = "variable", values_to = "value") %>%
  filter(variable != "total_pop", !str_starts(variable, "total_")) %>%
  mutate(
    day   = readr::parse_number(variable),
    type  = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
    weeks = ceiling((day-1) / 7) # -1 needed to correct day counter indexing [important for counting last week correctly, 546 total days is divisible by 7]
  ) %>%
  group_by(weeks, type, age_group) %>%
  summarise(cases = sum(value), .groups = "drop") %>%
  left_join(
    results %>%
      group_by(age_group) %>%
      summarise(total_pop = sum(total_pop)),
    by = "age_group"
  ) %>%
  mutate(inc = cases / total_pop * 100000) %>%
  dplyr::select(-total_pop)

#write.csv(age_group_results, paste0("stochasticity-results-012026-18months-CF65/",
#                                    "seed", job_selector, ".csv"))


###################################################################################################
# PLOTTING
# Read in observed data
observed_severe_incidence <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput-updated.csv")[,-1] %>%
  filter(weeks_since %in% c(0:95)) %>% # filter(weeks_since %in% c(0:79)) %>% # KMB: choose which one
  dplyr::select(age_group, weeks_since, adj_inc) %>% 
  rename(inc = adj_inc, weeks = weeks_since) %>%
  mutate(age_group = case_when(
    age_group == "≥75 years" ~ "75+ years",
    age_group == "0-17 years (Children)" ~ "0-17 years",
    TRUE ~ age_group
  ))

combined <- rbind(age_group_results %>% 
                    filter(type == "severe") %>% 
                    dplyr::select(-type, -cases) %>%
                    mutate(group = "Simulated"),
                    observed_severe_incidence %>% 
                    #filter(weeks_since %in% c(0:79)) %>% 
                    #rename( inc = adj_inc) %>% 
                    dplyr::select(age_group, weeks, inc) %>% 
                    mutate(group = "Observed"))

# Create plot with all 100 simulations
ggplot(data = combined, aes(x = weeks, y = inc, color = age_group, group = age_group)) + 
  geom_line(data = combined %>% filter(group == "Simulated"), alpha = 0.4, linewidth = 0.5) +
  geom_line(data = combined %>% filter(group == "Observed"), linewidth = 1) +
  geom_point(data = combined %>% filter(group == "Observed"), size = 2) +
  ylab("Severe Incidence (per 100,000 persons)") +
  xlab("Time (weeks)") +
  #ylim(0, 35) + 
  #xlim(0, 78) + 
  #ggtitle() +
  theme(title = element_text(size = 10))
