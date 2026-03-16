###################################################################################################
# Title: Run ABC-SMC with EasyABC
# Author: Kate Bubar
# Date: February 5, 2026
###################################################################################################

library(here)
library(data.table)
library(dplyr)
library(EasyABC)
library(magrittr)

# Get the number of cores from the SLURM environment variable
n_cores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE"))

# Print to confirm (will appear in your log file)
cat("Using", n_cores, "cores\n")

#===============================================================================
# CONFIGURATION
#===============================================================================
lambda_grid <- c(0.5, 0.4, 0.45, 0.45, 0.55, 0.65, 0.55, 0.6, 0.55, 
                 0.5, 0.5, 0.5,  0.55, 0.5,  0.6,  0.55, 0.6, 0.65)

# Number of cores for parallelization
# n_cores <- 3  # Set to 1 for serial, or detectCores() - 2 for local

# # Prior bounds
# prior_low <- c(0.425, 0.35, 0.4, 0.4, 0.5, 0.55, 0.5, 0.006)
# prior_upp <- c(0.55, 0.5, 0.5, 0.55, 0.6, 0.65, 0.65, 0.018)
prior_low <- c(lambda_grid[1:9]*0.85, 0.006)
prior_upp <- c(lambda_grid[1:9]*1.15, 0.018)

# ABC-SMC settings
N <- 50  # Number of particles per generation
sim_days <- 274
sim_months <- 9

#===============================================================================
# OBSERVED SUMMARY STATISTICS (target)
#===============================================================================
observed_severe <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput.csv")[,-1] %>%
  filter(weeks_since %in% c(0:79), 
         age_group != "0-17 years") %>%
  dplyr::select(age_group, weeks_since, adj_inc) %>% 
  rename(weeks = weeks_since, inc = adj_inc) %>% 
  arrange(weeks, age_group) %>%
  mutate(age_group = case_when(
    age_group == "≥75 years" ~ "75+ years",
    age_group == "0-17 years (Children)" ~ "0-17 years",
    TRUE ~ age_group
  )) %>%
  filter(weeks <= floor(sim_days/7))

# Calculate weekly totals for 75+
observed_weekly_75 <- observed_severe$inc[observed_severe$age_group == "75+ years"]
observed_weekly_65_74 <- observed_severe$inc[observed_severe$age_group == "65-74 years"]
observed_weekly_50_64 <- observed_severe$inc[observed_severe$age_group == "50-64 years"]
observed_weekly_30_49 <- observed_severe$inc[observed_severe$age_group == "30-49 years"]
observed_weekly_18_29 <- observed_severe$inc[observed_severe$age_group == "18-29 years"]

# Convert month boundaries (in days) to week indices
# Day 1 = Week 1, Day 7 = Week 1, Day 8 = Week 2, etc.
month_boundaries <- c(1, 32, 63, 93, 124, 154, 185, 216, 244, 275, 305, 336, 
                      366, 397, 428, 458, 489, 519, 548)
month_boundaries_weeks <- ceiling(month_boundaries / 7)

# Aggregate observed weekly data using same month boundaries
aggregate_to_monthly <- function(weekly_data) {
  sapply(1:sim_months, function(i) {
    week_start <- month_boundaries_weeks[i]
    week_end <- month_boundaries_weeks[i + 1] - 1
    weeks_in_month <- week_start:week_end
    weeks_in_month <- weeks_in_month[weeks_in_month >= 1 & weeks_in_month <= 79]
    sum(weekly_data[weeks_in_month])
  })
}

# Apply to all age groups
observed_monthly_18_29 <- aggregate_to_monthly(observed_weekly_18_29)
observed_monthly_30_49 <- aggregate_to_monthly(observed_weekly_30_49)
observed_monthly_50_64 <- aggregate_to_monthly(observed_weekly_50_64)
observed_monthly_65_74 <- aggregate_to_monthly(observed_weekly_65_74)
observed_monthly_75 <- aggregate_to_monthly(observed_weekly_75)

# Calculate cumulative incidence by age group
observed_cumul <- observed_severe %>%
  group_by(age_group) %>%
  summarize(cumulative_inc = sum(inc, na.rm = TRUE), .groups = 'drop')

summary_stat_obs <- c(
  observed_monthly_18_29,
  observed_monthly_30_49,
  observed_monthly_50_64,
  observed_monthly_65_74,
  observed_monthly_75
  #observed_weekly_75,
  
  # # Age-specific burden
  # cumulative_incidence_18_29 = observed_cumul$cumulative_inc[observed_cumul$age_group == "18-29 years"],
  # cumulative_incidence_30_49 = observed_cumul$cumulative_inc[observed_cumul$age_group == "30-49 years"],
  # cumulative_incidence_50_64 = observed_cumul$cumulative_inc[observed_cumul$age_group == "50-64 years"],
  # cumulative_incidence_65_74 = observed_cumul$cumulative_inc[observed_cumul$age_group == "65-74 years"],
  # cumulative_incidence_75 = observed_cumul$cumulative_inc[observed_cumul$age_group == "75+ years"]
  
  # # Temporal pattern (use 75+ as an example)
  # peak_incidence_75 = max(observed_weekly_75),
  # peak_week_75 = which.max(observed_weekly_75),
  # mean_incidence_75 = mean(observed_weekly_75)
)

#===============================================================================
# MODEL WRAPPER FOR EasyABC
#===============================================================================

# This function is what EasyABC will call for each parameter set
run_model_for_abc <- function(param) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  
  set.seed(param[1])
  
  # Set parameters
  lambda_1 <<- param[2]
  lambda_2 <<- param[3]
  lambda_3 <<- param[4]
  lambda_4 <<- param[5]
  lambda_5 <<- param[6]
  lambda_6 <<- param[7]
  lambda_7 <<- param[8]
  lambda_8 <<- param[9]
  lambda_9 <<- param[10]

  baseline_case_hosp_frac <<- param[11]
  
  # Set other lambdas to NA
  lambda_10 <<- NA #param[10]
  lambda_11 <<- NA #param[11]
  lambda_12 <<- NA #param[12]
  lambda_13 <<- NA #param[13]
  lambda_14 <<- NA #param[14]
  lambda_15 <<- NA #param[15]
  lambda_16 <<- NA #param[16]
  lambda_17 <<- NA #param[17]
  lambda_18 <<- NA #param[18]
  
  time_since <<- 0
  sim_days <<- 274
  
  # Source inputs (needs variables in environment)
  #source(here::here("functions/vaccine-uptake-assignment-historical.R"))
  source(here::here("functions/simulation-functions-gridsearch.R"))
  source(here::here("functions/simulation-inputs-historical-5mil-preallocatedvax.R"))
  
  # Run simulation
  set.seed(1)
  results1 <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = sim_days, rel_ve_check = FALSE)
  
  set.seed(2)
  results2 <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = sim_days, rel_ve_check = FALSE)
  
  set.seed(3)
  results3 <- simulation_semiannual_strat_8_9_10(df_sim, sim_days = sim_days, rel_ve_check = FALSE)
  
  # Function to process results
  process_results <- function(results) {
    age_group_results <- results %>%
      group_by(age_group) %>%
      summarise(across(-risk_group, sum)) %>%
      pivot_longer(-age_group, names_to = "variable", values_to = "value") %>%
      filter(variable != "total_pop", 
             !str_starts(variable, "total_")) %>%
      mutate(
        day   = readr::parse_number(variable),
        type  = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
        weeks = ceiling((day - 1) / 7)
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
      dplyr::select(-total_pop) %>%
      filter(age_group != "0-17 years", type == "severe", weeks <= floor(sim_days/7)) %>%
      dplyr::select(-cases, -type)
    
    list(
      cumul = age_group_results %>%
        group_by(age_group) %>%
        summarize(cumulative_inc = sum(inc, na.rm = TRUE), .groups = 'drop'),
      weekly = age_group_results %>%
        filter(age_group == "75+ years") %>%
        arrange(weeks)
    )
  }
  
  # Process and average
  processed <- lapply(list(results1, results2, results3), process_results)
  
  simulated_cumul <- bind_rows(lapply(processed, `[[`, "cumul")) %>%
    group_by(age_group) %>%
    summarize(cumulative_inc = mean(cumulative_inc), .groups = 'drop')
  
  simulated_weekly <- bind_rows(lapply(processed, `[[`, "weekly")) %>%
    group_by(weeks) %>%
    summarize(inc = mean(inc), .groups = 'drop')
  
  # age_group_results <- results %>%
  #   group_by(age_group) %>%
  #   summarise(across(-risk_group, sum)) %>%
  #   pivot_longer(-age_group, names_to = "variable", values_to = "value") %>%
  #   filter(variable != "total_pop", 
  #          !str_starts(variable, "total_")) %>%
  #   mutate(
  #     day   = readr::parse_number(variable),
  #     type  = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
  #     weeks = ceiling((day - 1) / 7)
  #   ) %>%
  #   group_by(weeks, type, age_group) %>%
  #   summarise(cases = sum(value), .groups = "drop") %>%
  #   left_join(
  #     results %>%
  #       group_by(age_group) %>%
  #       summarise(total_pop = sum(total_pop)),
  #     by = "age_group"
  #   ) %>%
  #   mutate(inc = cases / total_pop * 100000) %>%
  #   dplyr::select(-total_pop) %>%
  #   filter(age_group != "0-17 years") %>%
  #   filter(type == "severe") %>%
  #   dplyr::select(-cases, -type) %>%
  #   filter(weeks <= floor(sim_days/7))
  # 
  # simulated_cumul <- age_group_results %>%
  #   group_by(age_group) %>%
  #   summarize(cumulative_inc = sum(inc, na.rm = TRUE), .groups = 'drop')
  # 
  # simulated_weekly <- age_group_results %>%
  #   filter(age_group == "75+ years") %>%
  #   arrange(weeks)
  
  # # Convert month boundaries (in days) to week indices
  # # Day 1 = Week 1, Day 7 = Week 1, Day 8 = Week 2, etc.
  # month_boundaries <- c(1, 32, 63, 93, 124, 154, 185, 216, 244, 275, 305, 336, 
  #                       366, 397, 428, 458, 489, 519, 548)
  # month_boundaries_weeks <- ceiling(month_boundaries / 7)
  # 
  # # Aggregate observed weekly data using same month boundaries
  # simulated_monthly_75 <- sapply(1:18, function(i) {
  #   week_start <- month_boundaries_weeks[i]
  #   week_end <- month_boundaries_weeks[i + 1] - 1
  #   
  #   # Keep only weeks that exist in your data
  #   weeks_in_month <- week_start:week_end
  #   weeks_in_month <- weeks_in_month[weeks_in_month >= 1 & weeks_in_month <= 79]
  #   
  #   sum(simulated_weekly$inc[weeks_in_month])
  # })
  
  results <- c(
    simulated_weekly$inc,
    #simulated_monthly_75,
    
    # cumulative_incidence_18_29 = simulated_cumul$cumulative_inc[simulated_cumul$age_group == "18-29 years"],
    # cumulative_incidence_30_49 = simulated_cumul$cumulative_inc[simulated_cumul$age_group == "30-49 years"],
    # cumulative_incidence_50_64 = simulated_cumul$cumulative_inc[simulated_cumul$age_group == "50-64 years"],
    # cumulative_incidence_65_74 = simulated_cumul$cumulative_inc[simulated_cumul$age_group == "65-74 years"],
    # cumulative_incidence_75 = simulated_cumul$cumulative_inc[simulated_cumul$age_group == "75+ years"]
    # 
    # peak_incidence_75 = max(simulated_weekly$inc, na.rm = TRUE),
    # peak_week_75 = which.max(simulated_weekly$inc),
    # mean_incidence_75 = mean(simulated_weekly$inc, na.rm = TRUE)
  )
  
  return(results)
}

# start_time <- Sys.time()
# run_model_for_abc(par)
# end_time <- Sys.time()
# 
# cat(sprintf("Total time: %.1f minutes\n", 
#             as.numeric(end_time - start_time, units = "mins")))


#===============================================================================
# DEFINE PRIOR
#===============================================================================
# EasyABC prior format: list of sampling functions
prior_functions <- list(c("unif", prior_low[1], prior_upp[1]), #lambda1
                        c("unif", prior_low[2], prior_upp[2]), #lambda2
                        c("unif", prior_low[3], prior_upp[3]), #lambda3
                        c("unif", prior_low[4], prior_upp[4]), #lambda4
                        c("unif", prior_low[5], prior_upp[5]), #lambda5
                        c("unif", prior_low[6], prior_upp[6]), #lambda6
                        c("unif", prior_low[7], prior_upp[7]), #lambda7
                        c("unif", prior_low[8], prior_upp[8]), #lambda8
                        c("unif", prior_low[9], prior_upp[9]), #lambda9
                        # c("unif", prior_low[10], prior_upp[10]), #lambda10
                        # c("unif", prior_low[11], prior_upp[11]), #lambda11
                        # c("unif", prior_low[12], prior_upp[12]), #lambda12
                        # c("unif", prior_low[13], prior_upp[13]), #lambda13
                        # c("unif", prior_low[14], prior_upp[14]), #lambda14
                        # c("unif", prior_low[15], prior_upp[15]), #lambda15
                        # c("unif", prior_low[16], prior_upp[16]), #lambda16
                        # c("unif", prior_low[17], prior_upp[17]), #lambda17
                        # c("unif", prior_low[18], prior_upp[18]), #lambda18
                        c("unif", prior_low[10], prior_upp[10])  #baseline CHF
                        ) 


#===============================================================================
# RUN ABC-SMC
#===============================================================================

start_time <- Sys.time()

abc_results <- ABC_sequential(
  method = "Delmoral",  # Standard ABC-SMC algorithm
  model = run_model_for_abc,
  alpha = 0.9, # 0.9 is default for delmoral
  prior = prior_functions,
  #prior_test = "abs(X1 - X2) < 1", # && abs(X2 - X3) < 0.15 && abs(X3 - X4) < 0.15 && abs(X4 - X5) < 0.15 && abs(X5 - X6) < 0.15", #&& abs(X6 - X7) < 0.15 && abs(X7 - X8) < 0.15 && abs(X8 - X9) < 0.15 && abs(X9 - X10) < 0.15 && abs(X10 - X11) < 0.15 && abs(X11 - X12) < 0.15 && abs(X12 - X13) < 0.15 && abs(X13 - X14) < 0.15 && abs(X14 - X15) < 0.15 && abs(X15 - X16) < 0.15 && abs(X16 - X17) < 0.15 && abs(X17 - X18) < 0.15", 
  nb_simul = N, # initial number of accepted simulations of the model
  M = 1,
  summary_stat_target = summary_stat_obs,
  #dist_weights = c(rep(1, 5), 2, 2, 1),
  tolerance_target = 100, # used to be 10
  use_seed = TRUE, # was true
  n_cluster = n_cores,
  verbose = TRUE
)

# abc_results <- ABC_sequential(
#   method = "Lenormand",  # Standard ABC-SMC algorithm
#   model = run_model_for_abc,
#   alpha = 0.5, # 0.5 is default for Lenormand
#   prior = prior_functions,
#   #prior_test = "abs(X1 - X2) < 0.1 && abs(X2 - X3) < 0.1 && abs(X3 - X4) < 0.1 && abs(X4 - X5) < 0.1 && abs(X5 - X6) < 0.1",
#   nb_simul = N, # initial number of accepted simulations of the model
#   summary_stat_target = summary_stat_obs,
#   p_acc_min = 0.05, # stopping criterion
#   #use_seed = TRUE,
#   #n_cluster = n_cores,
#   verbose = TRUE
#   )

# && abs(X6 - X7) < 0.1 &&
#abs(X7 - X8) < 0.1 && abs(X8 - X9) < 0.1 && abs(X9 - X10) < 0.1 &&
#abs(X10 - X11) < 0.1 && abs(X11 - X12) < 0.1 && abs(X12 - X13) < 0.1 &&
#abs(X13 - X14) < 0.1 && abs(X14 - X15) < 0.1 && abs(X15 - X16) < 0.1 &&
#abs(X16 - X17) < 0.1 && abs(X17 - X18) < 0.1

# abc_results$param: matrix of accepted parameter values in final generation
# abc_results$stats: summary statistics for the accepted sim. (same dim as simulated summary statistics)
# abc_results$weights: importance weights for each accepted particle
# abc_results$stats_normalization: The standard deviation of the summary statistics across the model simulations of the initial step. These values are used to normalize the summary statistics before the computation of the Euclidean distance between simulations and data.
# abc_results$epsilon: the tolerance threshold(s) used
# abc_results$nsim: total number of simulations performed across all generations

end_time <- Sys.time()

saveRDS(abc_results, paste0("abc_results_022426.rds"))

cat(sprintf("Total time: %.1f minutes\n", 
            as.numeric(end_time - start_time, units = "mins")))



