########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: March 4th, 2025
########################################################################################################################

rm(list=ls())
gc()

#Load libraries
library(tidyverse)
library(tibble)
library(lubridate)
library(here)
library(data.table)
library(reshape2)
library(dtw)

# Observed data (weekly severe incidence estimates, by age group) -- clean this up
observed_data <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod.csv")[,-1] %>%
  filter(weeks_since %in% c(0:77), age_group != "0-17 years") %>%
  dplyr::select(weeks_since, age_group, adj_inc) %>%
  arrange(weeks_since, age_group) %>% group_by(age_group) %>%
  pivot_wider(names_from = age_group, values_from = adj_inc) 

#Set simulation size 
sim_size_max <- 10000

#Set up folder structure to save simulation results
dir.create("mcmc-output")

#First checking if output file exists to pick up mcmc where it left off
# NOTE: Make sure that the output file name is the same one in `mcmc-functions.R` script
if (file.exists("mcmc-output/mcmc-params-18lambdas-window2-052325.csv")) {
  
  #Read in output file and filter saved params, pick up mcmc from here
  output <- read.csv("mcmc-output/mcmc-params-18lambdas-window2-052325.csv")[,-1] %>% 
    filter(!is.na(V1))
  
  last_params <- as.numeric(output %>% tail(1))
  
  #Set initial params
  starting_params <- c(lambda_1 = last_params[1],
                       lambda_2 = last_params[2],
                       lambda_3 = last_params[3],
                       lambda_4 = last_params[4],
                       lambda_5 = last_params[5],
                       lambda_6 = last_params[6],
                       lambda_7 = last_params[7],
                       lambda_8 = last_params[8],
                       lambda_9 = last_params[9],
                       lambda_10 = last_params[10],
                       lambda_11 = last_params[11],
                       lambda_12 = last_params[12],
                       lambda_13 = last_params[13],
                       lambda_14 = last_params[14],
                       lambda_15 = last_params[15],
                       lambda_16 = last_params[16],
                       lambda_17 = last_params[17],
                       lambda_18 = last_params[18],
                       baseline_case_hosp_frac = last_params[19],
                       time_since = last_params[20])
  
  print(starting_params)

  #index of most recently saved params
  last_sim_index <- nrow(output)
  
  # Reading model scripts
  source(here::here("model-setup-onetime-1mil.R"))
  source(here::here("model-functions.R"))
  source(here::here("mcmc-functions.R"))
  
  
  #Run MCMC
  chain = run_metropolis_MCMC(starting_params, sim_size_max, last_sim_index)
  
    
} else {
  
  #Otherwise set initial params (best guess)
  starting_params <- c(lambda_1 = 1.00483905470604,
                       lambda_2 = 0.999239300669564,
                       lambda_3 = 1.00131677474849,
                       lambda_4 = 1.00171917687992,
                       lambda_5 = 1.00865700800111,
                       lambda_6 = 1.00029693936008,
                       lambda_7 = 0.997286977916913,
                       lambda_8 = 0.996986977916913,
                       lambda_9 = 1.00400330011814,
                       lambda_10 = 0.997251454980608,
                       lambda_11 = 1.001,
                       lambda_12 = 0.997204999183731,
                       lambda_13 = 1.00023367016819,
                       lambda_14 = 1.00331336257227,
                       lambda_15 = 0.997423744845859,
                       lambda_16 = 1.00302936138868 ,
                       lambda_17 = 1.001,
                       lambda_18 = 1.0016,
                       baseline_case_hosp_frac = 0.00542120442934245,
                       time_since = -7.4864)

  # Reading model scripts
  source(here::here("model-setup-onetime-1mil.R"))
  source(here::here("model-functions.R"))
  source(here::here("mcmc-functions.R"))

  
  #Run MCMC
  chain = run_metropolis_MCMC(starting_params, sim_size_max, 1)
  
  
}

#Create an empty file to signal to Marlowe that the job is complete
file.create("job-complete-18lambdas-window2-meanUI-052325.txt")

