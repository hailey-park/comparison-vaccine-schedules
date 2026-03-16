###################################################################################################
# Title: Run ABC-SMC with EasyABC
# Author: Kate Bubar
# Date: February 5, 2026
###################################################################################################
rm(list = ls())

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
                 0.5, 0.5, 0.5, 0.55, 0.5, 0.6, 0.55, 0.6, 0.65)

# Number of cores for parallelization
# n_cores <- 3  # Set to 1 for serial, or detectCores() - 2 for local

# # Prior bounds
# prior_low <- c(0.425, 0.35, 0.4, 0.4, 0.5, 0.55, 0.5, 0.006)
# prior_upp <- c(0.55, 0.5, 0.5, 0.55, 0.6, 0.65, 0.65, 0.018)
prior_low <- c(lambda_grid[1:7]*0.85, 0.006)
prior_upp <- c(lambda_grid[1:7]*1.15, 0.018)

# ABC-SMC settings
N <- 20  # Number of particles per generation
sim_days <- 215 # 7 months

# Output
#output_dir <- "output/calibration"
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

# Calculate cumulative incidence by age group
observed_cumul <- observed_severe %>%
  group_by(age_group) %>%
  summarize(cumulative_inc = sum(inc, na.rm = TRUE), .groups = 'drop')

# Calculate weekly totals for 75+
observed_weekly_75 <- observed_severe[observed_severe$age_group == "75+ years",]$inc

summary_stat_obs <- c(
  observed_weekly_75,
  
  # Age-specific burden
  cumulative_incidence_18_29 = observed_cumul[observed_cumul$age_group == "18-29 years",]$cumulative_inc,
  cumulative_incidence_30_49 = observed_cumul[observed_cumul$age_group == "30-49 years",]$cumulative_inc,
  cumulative_incidence_50_64 = observed_cumul[observed_cumul$age_group == "50-64 years",]$cumulative_inc,
  cumulative_incidence_65_74 = observed_cumul[observed_cumul$age_group == "65-74 years",]$cumulative_inc,
  cumulative_incidence_75 = observed_cumul[observed_cumul$age_group == "75+ years",]$cumulative_inc

  # # Temporal pattern (use 75+ as an example)
  # peak_incidence_75 = max(observed_weekly_75),
  # peak_week_75 = which.max(observed_weekly_75),
  # mean_incidence_75 = mean(observed_weekly_75)
)

#===============================================================================
# MODEL WRAPPER FOR EasyABC
#===============================================================================

# This function is what EasyABC will call for each parameter set
run_model_for_abc <- function(par) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  
  set.seed(par[1])
  
  # Set parameters
  lambda_1 <<- par[2]
  lambda_2 <<- par[3]
  lambda_3 <<- par[4]
  lambda_4 <<- par[5]
  lambda_5 <<- par[6]
  lambda_6 <<- par[7]
  lambda_7 <<- par[8]
  baseline_case_hosp_frac <<- par[9]
  
  # Set other lambdas to NA
  lambda_8 <<- lambda_9 <<- lambda_10 <<- lambda_11 <<- lambda_12 <<- NA
  lambda_13 <<- lambda_14 <<- lambda_15 <<- lambda_16 <<- lambda_17 <<- lambda_18 <<- NA
  
  time_since <<- 0
  sim_days <<- 215
  
  # Source inputs (needs variables in environment)
  source(here::here("functions/vaccine-uptake-assignment-historical.R"))
  source(here::here("functions/simulation-functions-gridsearch.R"))
  source(here::here("functions/simulation-inputs-historical-5mil.R"))
  
  # Run simulation
  df_sim <- historical_vax_assignment(clean_df[[1]], first_dose_coverage, 
                                      second_dose_coverage, next_year_dose_data)
  
  results <- simulation_semiannual_strat_8_9_10(df_sim, 
                                                sim_days = sim_days, 
                                                rel_ve_check = FALSE)
  
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
    filter(age_group != "0-17 years") %>%
    filter(type == "severe") %>%
    dplyr::select(-cases, -type) %>%
    filter(weeks <= floor(sim_days/7))
      
  simulated_cumul <- age_group_results %>%
    group_by(age_group) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE), .groups = 'drop')
  
  simulated_weekly <- age_group_results %>%
    filter(age_group == "75+ years") %>%
    arrange(weeks)
  
  results <- c(
    simulated_weekly$inc,
    cumulative_incidence_18_29 = simulated_cumul[simulated_cumul$age_group == "18-29 years",]$cumulative_inc,
    cumulative_incidence_30_49 = simulated_cumul[simulated_cumul$age_group == "30-49 years",]$cumulative_inc,
    cumulative_incidence_50_64 = simulated_cumul[simulated_cumul$age_group == "50-64 years",]$cumulative_inc,
    cumulative_incidence_65_74 = simulated_cumul[simulated_cumul$age_group == "65-74 years",]$cumulative_inc,
    cumulative_incidence_75 = simulated_cumul[simulated_cumul$age_group == "75+ years",]$cumulative_inc
    # 
    # peak_incidence_75 = max(simulated_weekly$inc, na.rm = TRUE),
    # peak_week_75 = which.max(simulated_weekly$inc),
    # mean_incidence_75 = mean(simulated_weekly$inc, na.rm = TRUE)
  )
  
  return(results)
}

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
                        c("unif", prior_low[8], prior_upp[8])) #baseline CHF


# prior_functions <- list(
#   c("normal", lambda_grid[1], 0.15),
#   c("normal", lambda_grid[2], 0.15),
#   c("normal", lambda_grid[3], 0.15),
#   c("normal", lambda_grid[4], 0.15),
#   c("normal", lambda_grid[5], 0.15),
#   c("normal", lambda_grid[6], 0.15),
#   c("unif", prior_low[7], prior_upp[7])
# )
start_time <- Sys.time()
run_model_for_abc(par)
end_time <- Sys.time()

cat(sprintf("Total time: %.1f minutes\n", 
            as.numeric(end_time - start_time, units = "mins")))



#===============================================================================
# RUN ABC-SMC
#===============================================================================

start_time <- Sys.time()

# # Run ABC-SMC
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

abc_results <- ABC_sequential(
  method = "Delmoral",  # Standard ABC-SMC algorithm
  model = run_model_for_abc,
  alpha = 0.9, # 0.9 is default for delmoral
  prior = prior_functions,
  prior_test = "abs(X1 - X2) < 0.1 && abs(X2 - X3) < 0.1 && abs(X3 - X4) < 0.1 && abs(X4 - X5) < 0.1 && abs(X5 - X6) < 0.1",
  nb_simul = N, # initial number of accepted simulations of the model
  summary_stat_target = summary_stat_obs,
  #dist_weights = c(rep(1, 5), 2, 2, 1),
  tolerance_target = 10, # used to be 10
  use_seed = TRUE,
  n_cluster = n_cores,
  verbose = TRUE
  )

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

saveRDS(abc_results, paste0("abc_results_job_021926.rds"))

cat(sprintf("Total time: %.1f minutes\n", 
            as.numeric(end_time - start_time, units = "mins")))



