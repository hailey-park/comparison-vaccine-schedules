###################################################################################################
#Title: Run Model - multiple iterations to understand stochasticity, designed to run locally
#Original author: Hailey Park
#Date: September 6, 2025
#Last Updated: January 26, 2026
###################################################################################################

rm(list=ls())

#start <- Sys.time()
#Load libraries
library(tidyverse)
library(here)
library(data.table)

###################################################################################################
#LOAD in scripts (only needs to be done once)
setwd("~/comparison-vaccine-schedules")
source("functions/vaccine-uptake-assignment-historical.R")
#source("functions/vaccine-uptake-assignment-scenarios.R") # FIXME - REMOVE FOR HISTORICAL
source("functions/simulation-functions-gridsearch.R") # includes how long to run model
#source("functions/simulation-functions-updated.R") # FIXME - REMOVE FOR HISTORICAL

#SETUP folder structure to save simulation results
folder <- "output/stochasticity-results-031626"

if (!dir.exists(folder)) {
  dir.create(folder)
}

###################################################################################################
#RUN a sherlock job of each parameter combination

#parameter_combos <- read.csv("sample_of_accepted_parameters_030626.csv")

total_runs <- c(8:10) # originally up to 10

#Get parameters specific to this_run
for (this_run in total_runs) {
  print(this_run)
  
  # Burn-in period: days 1-304
  sim_days <- 123
    
  lambda_1 <- 0.49
  lambda_2 <- 0.43
  lambda_3 <- 0.44
  lambda_4 <- 0.46
  lambda_5 <- 0.59
  lambda_6 <- 0.66
  lambda_7 <- 0.60
  lambda_8 <- 0.60
  lambda_9 <- 0.57
  lambda_10 <- 0.57
  
  # Phase 1
  lambda_11 <- 0.52-0.05
  lambda_12 <- 0.52-0.05
  lambda_13 <- 0.52-0.05 
  lambda_14 <- 0.52-0.05
  lambda_15 <- 0.52-0.05
  lambda_16 <- 0.55-0.05
  
  # Phase 2
  lambda_17 <- 0.6-0.05
  lambda_18 <- 0.65-0.05
  lambda_19 <- 0.5
  lambda_20 <- 0.5
  lambda_21 <- 0.5
  lambda_22 <- 0.5
  
  baseline_case_hosp_frac <- 0.0115
  
  time_since <- 0
  
  #RUN the simulation
  #   NOTE: vax distribution used in semiannual_strat_8_9_10 corresponds to the observed historical distributions
  
  # source("functions/simulation-inputs-historical-5mil-preallocatedvax.R") # FIXME - REMOVE FOR HISTORICAL
  
  start.time <- Sys.time()
  set.seed(this_run)
  source("functions/simulation-inputs-historical-5mil.R")
  df_sim <- historical_vax_assignment(clean_df[[1]], dose1_coverage, dose2_coverage, dose3_coverage, dose4_coverage)
  results <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = sim_days, rel_ve_check = FALSE)
  
  # # FIXME - REMOVE FOR HISTORICAL
  # source("functions/simulation-inputs-historical-5mil.R")
  # realistic_ind <- 1
  # df_sim <- strategy_6(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind)
  # results <- simulation_annual(df_sim)#, sim_days = sim_days, rel_ve_check = FALSE)

  print(paste0("job num: ", this_run, ", time: ", round(Sys.time() - start.time, 2)))
  
  age_group_summarised <- results %>%
    group_by(age_group) %>%
    summarise(across(c(-"risk_group"), sum))
  
  age_group_results <- age_group_summarised %>%
    pivot_longer(
      cols = -age_group,
      names_to = "variable",
      values_to = "cases"
    ) %>%
    filter(variable != "total_pop") %>%
    mutate(
      day   = readr::parse_number(variable),
      type  = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
      weeks = ceiling(day / 7)
    ) %>%
    group_by(weeks, type, age_group) %>%
    summarise(cases = sum(cases), .groups = "drop") %>%
    left_join(
      age_group_summarised %>%
        pivot_longer(
          cols = -age_group,
          names_to = "variable",
          values_to = "total_pop"
        ) %>%
         filter(variable == "total_pop") %>%
         dplyr::select(-variable),
      by = "age_group"
    ) %>%
    mutate(
      inc = cases / total_pop * 100000
    ) %>%
    dplyr::select(-total_pop)
  
  # store results for scenario analysis
  write.csv(results, paste0(folder, "/",
                                      "results", this_run, ".csv"))
  
  # store age_group_results for plotting
  write.csv(age_group_results, paste0(folder, "/",
                                      "age_group_results", this_run, ".csv"))
  
}

