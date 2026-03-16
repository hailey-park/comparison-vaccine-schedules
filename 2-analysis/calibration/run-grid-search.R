###################################################################################################
#Title: Run Grid Search
#Author: Hailey Park
#Date created: September 6, 2025
#Last updated: January 22, 2026 - KMB
###################################################################################################

rm(list=ls())
gc()

#start <- Sys.time()
#Load libraries
library(data.table)
library(tidyverse)
library(here)

args <- commandArgs(trailingOnly = TRUE)
job_selector <- (as.numeric(args[1])-1)
###################################################################################################

#Load in scripts (only needs to be done once)
source(here::here("functions/vaccine-uptake-assignment-calibration.R"))
source(here::here("functions/simulation-functions-updated-grid-search.R")) # UPDATE length of sim in this file
##source("simulation-functions-updated-t0.R")

# SETUP folder structure to save results
folder_path <- "grid-search-results-120725-18months" # UPDATE for current iteration

if (!dir.exists(paste0("output/calibration/", folder_path))) {
  dir.create(paste0("output/calibration/", folder_path))
}

###################################################################################################
#Run a Sherlock job of each parameter combination

# #Parameter combinations
job_nums <- 10

parameter_combos <- merge(data.frame(lambda_17 = seq(0.35, 1.0, by = 0.05)),
                          data.frame(lambda_18 = seq(0.35, 1.0, by = 0.05)), all = TRUE
) %>%
  mutate(job_num = 1:job_nums)


#Get parameters specific to job_selector
for (job_selector in c(2:job_nums)) { # for grid search
  
  specific_params <- parameter_combos %>% filter(job_num == job_selector)
  
  # lambda_1 <- specific_params$lambda_1
  # lambda_2 <- specific_params$lambda_2
  # lambda_3 <- specific_params$lambda_3
  # lambda_4 <- specific_params$lambda_4
  # lambda_5 <- specific_params$lambda_5
  # lambda_6 <- specific_params$lambda_6
  # lambda_7 <- specific_params$lambda_7
  lambda_17 <- specific_params$lambda_17
  lambda_18 <- specific_params$lambda_18
  
  baseline_case_hosp_frac <- specific_params$baseline_case_hosp_frac
  
  time_since <- 0
  
  #set.seed(942)
  set.seed(job_selector)
  
  if(!file.exists(paste0("output/calibration/", folder_path, "/",
                  "lambda_17", lambda_17,
                  "lambda_18", lambda_18, ".csv"))) {
  
  #Read in simulation inputs script and calculate inputs with new parameter set
  #source(here::here("simulation-inputs-historical-updated-1mil.R")) # FIXME KMB - add if statement about population size
  source((here::here("functions/simulation-inputs-historical-5mil.R")))
  
  # ASSIGN vaccine uptake and timing to individuals
  start.time <- Sys.time()
  df_sim <- historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind)
  print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
  
  # RUN simulation
  start.time <- Sys.time()
  results <- simulation_semiannual_strat_8_9_10(df = df_sim, sim_days = 547) # max sim_days is 547 days
  print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
  
  # SUMMARIZE results and save
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
        select(-variable),
      by = "age_group"
    ) %>%
    mutate(
      inc = cases / total_pop * 100000
    ) %>%
    select(-total_pop)
  
  #Write to csv
  write.csv(age_group_results, paste0("output/calibration/", folder_path, "/",
                                      #"lambda_1", lambda_1,
                                      # "lambda_2", lambda_2,
                                      # "lambda_3", lambda_3,
                                      "lambda_17", lambda_17,
                                      "lambda_18", lambda_18, ".csv"))
  
  print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
  
}

