###################################################################################################
#Title: Run Model - multiple iterations to understand stochasticity, designed to run on Sherlock
#Original author: Hailey Park
#Date: September 6, 2025
#Last Updated: January 26, 2026
###################################################################################################

rm(list=ls())
gc()

#start <- Sys.time()
#Load libraries
library(tidyverse)
library(here)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
job_selector <- (as.numeric(args[1])-1)
###################################################################################################

#LOAD in scripts (only needs to be done once)
source(here::here("functions/vaccine-uptake-assignment-historical.R"))
source(here::here("functions/simulation-functions-gridsearch.R")) # includes how long to run model

###################################################################################################
#RUN a sherlock job of each parameter combination

# First check if output already exists
output_file <- paste0("stochasticity-results/seed", job_selector, ".csv")

if (file.exists(output_file)) {
  cat(paste("Output file already exists for job", job_selector, "- skipping\n"))
} else {
  # Else load in the parameters
  parameter_combos <- read.csv("sample_of_accepted_parameters_022126.csv")
  job_nums <- length(parameter_combos)
  
  #Get parameters specific to job_selector
  if (job_selector %in% c(1:job_nums)) { # for grid search
    #for (job_selector in c(0:job_nums)) {
    print(job_selector)
    
    #start.time <- Sys.time()
    specific_params <- parameter_combos[(job_selector+1),]
    #if (job_selector < 10) {
    
    # lambda_1 <- 0.5
    # lambda_2 <- 0.4
    # lambda_3 <- 0.45
    # lambda_4 <- 0.45
    # lambda_5 <- 0.55
    # lambda_6 <- 0.65
    # lambda_7 <- 0.55
    # lambda_8 <- 0.6
    # lambda_9 <- 0.55
    # lambda_10 <- 0.5
    # lambda_11 <- 0.5
    # lambda_12 <- 0.5
    # lambda_13 <- 0.55
    # lambda_14 <- 0.5
    # lambda_15 <- 0.6
    # lambda_16 <- 0.55
    # lambda_17 <- 0.6
    # lambda_18 <- 0.65
    # baseline_case_hosp_frac <- 0.012 # 0.012
    
    lambda_1 <- specific_params$param_1
    lambda_2 <- specific_params$param_2
    lambda_3 <- specific_params$param_3
    lambda_4 <- specific_params$param_4
    lambda_5 <- specific_params$param_5
    lambda_6 <- specific_params$param_6
    lambda_7 <- specific_params$param_7
    lambda_8 <- 0.6
    lambda_9 <- 0.55
    lambda_10 <- 0.5
    lambda_11 <- 0.5
    lambda_12 <- 0.5
    lambda_13 <- 0.55
    lambda_14 <- 0.5
    lambda_15 <- 0.6
    lambda_16 <- 0.55
    lambda_17 <- 0.6
    lambda_18 <- 0.65
    baseline_case_hosp_frac <- specific_params$param_8
    
    time_since <- 0
    set.seed(job_selector)
    
    #RUN the simulation
    #   NOTE: vax distribution used in semiannual_strat_8_9_10 corresponds to the observed historical distributions
    
    start.time <- Sys.time()
    source(here::here("functions/simulation-inputs-historical-5mil.R"))
    #print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
    
    #start.time <- Sys.time()
    df_sim <- historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data)
    #print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
    
    #start.time <- Sys.time()
    # Run simulation
    results <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = 215, rel_ve_check = FALSE)
    print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
    
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
    
    write.csv(age_group_results, paste0("stochasticity-results/",
                                        "seed", job_selector, ".csv"))
    #print(paste0("job num: ", job_selector, ", time: ", round(Sys.time() - start.time, 2)))
    
  }
}

###################################################################################################