###################################################################################################
# Title: Run ABC-SMC (Simplified Parallel Version)
# Author: Kate Bubar
# Date: February 5, 2026
###################################################################################################

library(tidyverse)
library(here)
library(data.table)
library(tmvtnorm)
library(foreach)
library(doParallel)

#===============================================================================
# CONFIGURATION
#===============================================================================

# Detect available cores and use conservatively
max_cores <- parallel::detectCores()
n_cores <- max(1, max_cores - 11)  # Leave 2 cores free for system

# Load scripts
source(here::here("functions/vaccine-uptake-assignment-historical.R"))
source(here::here("functions/simulation-functions-gridsearch.R"))

# Prior bounds
prior_low <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.006)
prior_upp <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.018)

# ABC-SMC settings
N <- 2
n_par <- 7

# Tolerance
epsilon_weekly <- 9000
epsilon_cumul <- 43000

# Output
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===============================================================================
# DATA
#===============================================================================

observed_severe <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput.csv")[,-1] %>%
  filter(weeks_since %in% c(0:79), age_group != "0-17 years") %>%
  dplyr::select(age_group, weeks_since, adj_inc) %>% 
  rename(weeks = weeks_since, inc = adj_inc) %>% 
  arrange(weeks, age_group) %>%
  mutate(age_group = case_when(
    age_group == "≥75 years" ~ "75+ years",
    age_group == "0-17 years (Children)" ~ "0-17 years",
    TRUE ~ age_group
  ))

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

sample_prior <- function() {
  c(
    lambda_1 = runif(1, prior_low[1], prior_upp[1]),
    lambda_2 = runif(1, prior_low[2], prior_upp[2]),
    lambda_3 = runif(1, prior_low[3], prior_upp[3]),
    lambda_4 = runif(1, prior_low[4], prior_upp[4]),
    lambda_5 = runif(1, prior_low[5], prior_upp[5]),
    lambda_6 = runif(1, prior_low[6], prior_upp[6]),
    baseline_case_hosp_frac = runif(1, prior_low[7], prior_upp[7])
  )
}

prior_non_zero <- function(par) {
  all(par >= prior_low & par <= prior_upp)
}

#===============================================================================
# DISTANCE FUNCTION
#===============================================================================

calc_distance <- function(simulated_data, sim_days) {
  
  simulated_filtered <- simulated_data %>% 
    filter(age_group != "0-17 years") %>%
    filter(type == "severe") %>%
    select(-cases, -type)
  
  observed_filtered <- observed_severe %>%
    filter(weeks <= floor(sim_days/7))
  
  combined <- simulated_filtered %>%
    full_join(observed_filtered, 
              by = c("age_group", "weeks"),
              suffix = c("_sim", "_obs"))
  
  weekly_sse <- combined %>%
    summarise(sse = sum((inc_sim - inc_obs)^2, na.rm = TRUE)) %>%
    pull(sse)
  
  cumulative_sse <- combined %>%
    group_by(age_group) %>%
    summarise(
      cumul_sim = sum(inc_sim, na.rm = TRUE),
      cumul_obs = sum(inc_obs, na.rm = TRUE)
    ) %>%
    summarise(sse = sum((cumul_sim - cumul_obs)^2, na.rm = TRUE)) %>%
    pull(sse)
  
  return(c(weekly_sse, cumulative_sse))
}

#===============================================================================
# SINGLE SIMULATION WRAPPER
#===============================================================================

run_one_sim <- function(par) {
  # Set all parameters as global variables (needed for sourcing)
  lambda_1 <<- par[1]
  lambda_2 <<- par[2]
  lambda_3 <<- par[3]
  lambda_4 <<- par[4]
  lambda_5 <<- par[5]
  lambda_6 <<- par[6]
  baseline_case_hosp_frac <<- par[7]
  
  lambda_7 <<- lambda_8 <<- lambda_9 <<- lambda_10 <<- lambda_11 <<- lambda_12 <<- NA
  lambda_13 <<- lambda_14 <<- lambda_15 <<- lambda_16 <<- lambda_17 <<- lambda_18 <<- NA
  
  time_since <<- 0
  sim_days <<- 184
  
  # Source inputs
  source(here::here("functions/simulation-inputs-historical-5mil.R"))
  
  # Run simulation
  df_sim <- historical_vax_assignment(clean_df[[1]], first_dose_coverage, 
                                      second_dose_coverage, next_year_dose_data)
  
  results <- simulation_semiannual_strat_8_9_10(df_sim, 
                                                sim_days = sim_days, 
                                                rel_ve_check = FALSE)
  
  # Process results
  age_group_results <- results %>%
    group_by(age_group) %>%
    summarise(across(-risk_group, sum)) %>%
    pivot_longer(-age_group, names_to = "variable", values_to = "value") %>%
    filter(variable != "total_pop", !str_starts(variable, "total_")) %>%
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
    select(-total_pop)
  
  # Calculate distance
  dist <- calc_distance(age_group_results, sim_days)
  
  return(list(
    par = par,
    dist_weekly = dist[1],
    dist_cumul = dist[2],
    accepted = (dist[1] <= epsilon_weekly && dist[2] <= epsilon_cumul)
  ))
}

#===============================================================================
# MAIN ABC LOOP (SIMPLIFIED)
#===============================================================================

cat("\n=== Starting ABC-SMC Generation 1 ===\n")
cat(sprintf("Target: %d particles\n", N))
cat(sprintf("Cores: %d\n", n_cores))
cat(sprintf("Tolerance: weekly=%.0f, cumul=%.0f\n\n", epsilon_weekly, epsilon_cumul))

# Setup cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export everything to workers
clusterExport(cl, ls(), envir = environment())
clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)

# Load packages and source files on workers
clusterEvalQ(cl, {
  library(tidyverse)
  library(here)
  library(data.table)
  source(here::here("functions/vaccine-uptake-assignment-historical.R"))
  source(here::here("functions/simulation-functions-gridsearch.R"))
})

# Storage
accepted_particles <- list()
accepted_distances <- list()
total_sims <- 0

# Keep generating until we have N particles
while (length(accepted_particles) < N) {
  
  remaining <- N - length(accepted_particles)
  batch_size <- min(remaining, 500)  # Try 10x oversampling, max 500 at a time
  
  cat(sprintf("\nRunning batch of %d simulations (have %d/%d accepted)...\n", 
              batch_size, length(accepted_particles), N))
  
  start_time <- Sys.time()
  
  # Run simulations in parallel
  batch_results <- foreach(
    i = 1:batch_size,
    .packages = c('tidyverse', 'here', 'data.table'),
    .errorhandling = 'remove'  # Skip errors
  ) %dopar% {
    
    # Sample parameters
    par <- sample_prior()
    
    # Check bounds
    if (!prior_non_zero(par)) {
      return(NULL)
    }
    
    # Run simulation
    tryCatch({
      run_one_sim(par)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  end_time <- Sys.time()
  
  # Remove NULLs
  batch_results <- batch_results[!sapply(batch_results, is.null)]
  
  # Extract accepted
  newly_accepted <- batch_results[sapply(batch_results, function(x) x$accepted)]
  
  if (length(newly_accepted) > 0) {
    accepted_particles <- c(accepted_particles, newly_accepted)
  }
  
  # Update stats
  total_sims <- total_sims + batch_size
  n_accepted <- length(newly_accepted)
  acceptance_rate <- n_accepted / batch_size * 100
  
  cat(sprintf("  Completed in %.1f min\n", 
              as.numeric(end_time - start_time, units = "mins")))
  cat(sprintf("  Accepted: %d/%d (%.1f%%)\n", n_accepted, batch_size, acceptance_rate))
  cat(sprintf("  Total: %d/%d particles (%.1f%% overall)\n", 
              length(accepted_particles), N, 
              length(accepted_particles)/total_sims * 100))
  
  # Show distance stats if we have some
  if (length(newly_accepted) > 0) {
    weekly_dists <- sapply(newly_accepted, function(x) x$dist_weekly)
    cumul_dists <- sapply(newly_accepted, function(x) x$dist_cumul)
    cat(sprintf("  Weekly SSE: median=%.0f, range=[%.0f, %.0f]\n",
                median(weekly_dists), min(weekly_dists), max(weekly_dists)))
    cat(sprintf("  Cumul SSE:  median=%.0f, range=[%.0f, %.0f]\n",
                median(cumul_dists), min(cumul_dists), max(cumul_dists)))
  }
}

# Stop cluster
stopCluster(cl)

# Take first N particles
accepted_particles <- accepted_particles[1:N]

# Extract results
results_matrix <- do.call(rbind, lapply(accepted_particles, function(x) x$par))
dist_weekly <- sapply(accepted_particles, function(x) x$dist_weekly)
dist_cumul <- sapply(accepted_particles, function(x) x$dist_cumul)
weights <- rep(1/N, N)  # Equal weights for generation 1

# Save
results_df <- data.frame(
  lambda_1 = results_matrix[,1],
  lambda_2 = results_matrix[,2],
  lambda_3 = results_matrix[,3],
  lambda_4 = results_matrix[,4],
  lambda_5 = results_matrix[,5],
  lambda_6 = results_matrix[,6],
  baseline_case_hosp_frac = results_matrix[,7],
  weight = weights,
  dist_weekly = dist_weekly,
  dist_cumulative = dist_cumul
)

output_file <- file.path(output_dir, "loose_gen_1.csv")
write.csv(results_df, file = output_file, row.names = FALSE)

cat("\n=== Complete ===\n")
cat(sprintf("Total simulations: %d\n", total_sims))
cat(sprintf("Acceptance rate: %.2f%%\n", N/total_sims * 100))
cat(sprintf("Results saved to: %s\n", output_file))

# Print percentiles for next generation
cat("\n=== Percentiles for Gen 2 Tolerance ===\n")
cat(sprintf("Weekly:     25%%=%.0f, 50%%=%.0f, 75%%=%.0f\n",
            quantile(dist_weekly, 0.25),
            quantile(dist_weekly, 0.5),
            quantile(dist_weekly, 0.75)))
cat(sprintf("Cumulative: 25%%=%.0f, 50%%=%.0f, 75%%=%.0f\n",
            quantile(dist_cumul, 0.25),
            quantile(dist_cumul, 0.5),
            quantile(dist_cumul, 0.75)))
