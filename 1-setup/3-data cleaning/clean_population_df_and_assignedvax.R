# Clean the population dataframe (1 time)
# This will be read in in the simulation-inputs files
# Since we're not varying time_since anymore this is OK

time_since <- 0

#MAKE SURE YOU ARE READING IN THE CORRECT CALIBRATION FILE
entire_pop <- read.csv(paste0("data/clean-data/entire_population_model_initialization_daily_5mil_updated.csv"))[,-1] 

#clean population df
clean_population_df <- function(df){
  df %>% dplyr::select(c("individual", "age_group", "prior_inf", "risk_group","time_since_last_dose_inf", "time_since_last_dose", "prior_vacc")) %>%
    mutate(weeks_since_last_dose_inf = as.numeric(as.character(as.factor(interval(time_since_last_dose_inf, as.Date('2023-06-30')) %/% weeks(1)))) + 1,
           days_since_last_dose_inf = as.numeric(as.character(as.factor(interval(time_since_last_dose_inf, as.Date('2023-06-30')) %/% days(1)))),
           days_since_last_dose = as.numeric(as.character(as.factor(interval(time_since_last_dose, as.Date('2023-06-30')) %/% days(1)))),
           immuno = if_else(risk_group == "immunocompromised", 1, 0)) %>%
    rowwise() %>%
    mutate(weeks_since_last_dose_inf = max(min(weeks_since_last_dose_inf - ceiling(time_since), 104), 0),
           days_since_last_dose_inf = max(min(days_since_last_dose_inf - (ceiling(time_since) * 7), 730), 0),
           days_since_last_dose = min(days_since_last_dose, 730))
  
}

clean_df <- list(entire_pop) %>%
  lapply(clean_population_df)  

saveRDS(clean_df, "data/clean-data/clean_population_df_5mil.RDS")

rm(entire_pop)

# Save simulated pop with vaccines assigned for 5 mil population
first_dose_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/realistic-vaccine-coverage-by-age-and-risk.csv")[,-1] 
second_dose_coverage <- read.csv("data/clean-data/vax-uptake-scenarios/realistic-2nd-dose-annually-vaccine-coverage-by-age-and-risk.csv")[,-1]
next_year_dose_data <- read.csv("data/clean-data/vaccine-coverage-by-age-and-risk-forModelValidation.csv")[,-1]

start_time <- Sys.time()

df_sim <- historical_vax_assignment(clean_df[[1]], first_dose_coverage, 
                                    second_dose_coverage, next_year_dose_data)

end_time <- Sys.time()
cat(sprintf("Total time: %.1f minutes\n", 
            as.numeric(end_time - start_time, units = "mins")))

saveRDS(df_sim, "data/clean-data/df_sim_5mil.RDS")
