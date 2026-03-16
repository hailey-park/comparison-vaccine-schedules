########################################################################################################################
#Title: Main Analysis: Compare vaccination strategies
#Author: Hailey Park
#Date: July 1, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(reshape2)
library(data.table)
library(here)

args <- commandArgs(trailingOnly = TRUE)
job_selector <- (as.numeric(args[1])-1)

# #Calibrated baseline case-hospitalization fraction and time-since shift (load in)
# baseline_case_hosp_frac <- 0.005
# time_since <- --7
# 
# #Calibrated lambda parameters (load in)
# lambda_1 <- 1
# lambda_2 <- 0.999239300669564
# lambda_3 <- 1.00131677474849
# lambda_4 <- 1.00171917687992
# lambda_5 <- 1.00865700800111
# lambda_6 <- 1.00029693936008
# lambda_7 <- 0.997286977916913
# lambda_8 <- 0.996986977916913
# lambda_9 <- 1.00400330011814
# lambda_10 <- 0.997251454980608
# lambda_11 <- 1.001
# lambda_12 <- 0.997204999183731
# lambda_13 <- 1.00023367016819
# lambda_14 <- 1.00331336257227
# lambda_15 <- 0.997423744845859
# lambda_16 <- 1.00302936138868 
# lambda_17 <- 1.001
# lambda_18 <- 1.0016

# #Calibrated baseline case-hospitalization fraction and time-since shift (load in)
# baseline_case_hosp_frac <- 0.0105
# time_since <- 0

#Calibrated lambda parameters (load in)
lambda_1 <- 0.5
lambda_2 <- 0.4
lambda_3 <- 0.45
lambda_4 <- 0.45
lambda_5 <- 0.55
lambda_6 <- 0.65
lambda_7 <- 0.55
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
baseline_case_hosp_frac <- 0.012
time_since <- 0

#Call scripts
# NOTE: You might need to fix file paths for the script in order for the `here` function call to work
#setwd(here::here())
setwd("~/comparison-vaccine-schedules/analysis/4-main analysis")

print(getwd())
print(list.files())

# source(here::here("vaccine-uptake-assignment-scenarios.R"))
# source(here::here("simulation-functions-updated.R"))
source(("vaccine-uptake-assignment-scenarios.R"))
source(("simulation-functions-updated.R"))

#Set vaccine uptake assumption (realistic, optimistic)
if ((job_selector %in% c(1:120))) {
  #source(here::here("simulation-inputs-historical.R"))
  source(("simulation-inputs-historical.R"))
}

if ((job_selector %in% c(121:220))) {
  source(here::here("simulation-inputs-optimistic.R"))
}


#Create directory to save results
if (!dir.exists("final-analysis/simulation-results")) {
  dir.create("final-analysis/simulation-results", recursive = TRUE)
}

for (pop_strat in c("realistic", "optimistic")) {
  dir.create(paste0("final-analysis/simulation-results/", pop_strat))
}


if(job_selector %in% c(1:10)) {
  set.seed(job_selector %% 10)
  inspection_real <- simulation_semiannual_strat_8_9_10(historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_real, paste0("final-analysis/simulation-results/strat_real-updated-", job_selector %% 10, ".csv"))
}


if (job_selector %in% c(11:20)) {
  set.seed(job_selector %% 10)
  inspection_0 <- simulation_annual(strategy_0(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_0, paste0("final-analysis/simulation-results/realistic/strat_0-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(21:30)) {
  set.seed(job_selector %% 10)
  inspection_1 <- simulation_annual(strategy_1(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_1, paste0("final-analysis/simulation-results/realistic/strat_1-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(31:40)) {
  set.seed(job_selector %% 10)
  inspection_2 <- simulation_annual(strategy_2(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_2, paste0("final-analysis/simulation-results/realistic/strat_2-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(41:50)) {
  set.seed(job_selector %% 10)
  inspection_3 <- simulation_annual(strategy_3(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_3, paste0("final-analysis/simulation-results/realistic/strat_3-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(51:60)) {
  set.seed(job_selector %% 10)
  inspection_4 <- simulation_annual(strategy_4(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_4, paste0("final-analysis/simulation-results/realistic/strat_4-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(61:70)) {
  set.seed(job_selector %% 10)
  inspection_5 <- simulation_semiannual_strat_5_7(strategy_5(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_5, paste0("final-analysis/simulation-results/realistic/strat_5-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(71:80)) {
  set.seed(job_selector %% 10)
  inspection_6 <- simulation_annual(strategy_6(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_6, paste0("final-analysis/simulation-results/realistic/strat_6-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(81:90)) {
  set.seed(job_selector %% 10)
  inspection_7 <- simulation_semiannual_strat_5_7(strategy_7(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_7, paste0("final-analysis/simulation-results/realistic/strat_7-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(91:100)) {
  set.seed(job_selector %% 10)
  inspection_8 <- simulation_semiannual_strat_8_9_10(strategy_8(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_8, paste0("final-analysis/simulation-results/realistic/strat_8-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(101:110)) {
  set.seed(job_selector %% 10)
  inspection_9 <- simulation_semiannual_strat_8_9_10(strategy_9(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_9, paste0("final-analysis/simulation-results/realistic/strat_9-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(111:120)) {
  set.seed(job_selector %% 10)
  inspection_10 <- simulation_semiannual_strat_8_9_10(strategy_10(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_10, paste0("final-analysis/simulation-results/realistic/strat_10-", job_selector %% 10, ".csv"))
}


if (job_selector %in% c(121:130)) {
  set.seed(job_selector %% 10)
  inspection_1 <- simulation_annual(strategy_1(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_1, paste0("final-analysis/simulation-results/optimistic/strat_1-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(131:140)) {
  set.seed(job_selector %% 10)
  inspection_2 <- simulation_annual(strategy_2(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_2, paste0("final-analysis/simulation-results/optimistic/strat_2-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(141:150)) {
  set.seed(job_selector %% 10)
  inspection_3 <- simulation_annual(strategy_3(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_3, paste0("final-analysis/simulation-results/optimistic/strat_3-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(151:160)) {
  set.seed(job_selector %% 10)
  inspection_4 <- simulation_annual(strategy_4(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_4, paste0("final-analysis/simulation-results/optimistic/strat_4-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(161:170)) {
  set.seed(job_selector %% 10)
  inspection_5 <- simulation_semiannual_strat_5_7(strategy_5(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_5, paste0("final-analysis/simulation-results/optimistic/strat_5-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(171:180)) {
  set.seed(job_selector %% 10)
  inspection_6 <- simulation_annual(strategy_6(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_6, paste0("final-analysis/simulation-results/optimistic/strat_6-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(181:190)) {
  set.seed(job_selector %% 10)
  inspection_7 <- simulation_semiannual_strat_5_7(strategy_7(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_7, paste0("final-analysis/simulation-results/optimistic/strat_7-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(191:200)) {
  set.seed(job_selector %% 10)
  inspection_8 <- simulation_semiannual_strat_8_9_10(strategy_8(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_8, paste0("final-analysis/simulation-results/optimistic/strat_8-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(201:210)) {
  set.seed(job_selector %% 10)
  inspection_9 <- simulation_semiannual_strat_8_9_10(strategy_9(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_9, paste0("final-analysis/simulation-results/optimistic/strat_9-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(211:220)) {
  set.seed(job_selector %% 10)
  inspection_10 <- simulation_semiannual_strat_8_9_10(strategy_10(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_10, paste0("final-analysis/simulation-results/optimistic/strat_10-", job_selector %% 10, ".csv"))
}

