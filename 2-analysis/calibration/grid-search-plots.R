########################################################################################################################
#Title: Plotting Grid Search Runs
#Author: Hailey Park
#Date: January 16th, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(data.table)

setwd("/Users/katebubar/comparison-vaccine-schedules/analysis/3-model calibration/grid-search")
observed_severe_incidence_overall <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput.csv")[,-1] %>%
  filter(weeks_since %in% c(1:53)) %>% 
  dplyr::select(age_group, weeks_since, adj_inc) %>% 
  rename(observed_inc = adj_inc) %>%
  mutate(age_group = case_when(
    age_group == "≥75 years" ~ "75+ years",
    age_group == "0-17 years (Children)" ~ "0-17 years",
    TRUE ~ age_group
    ))

#setwd("/Users/katebubar/Desktop/gridsearch_results/grid-search-results-091925")
setwd("/Users/katebubar/Desktop/gridsearch_results/grid-search-results-092425-2months/grid-search-results-092425")

pdf("grid_search_plots_092425.pdf")

# Parameters that go with 091925
# for (lambda_1 in c(seq(0.9, 1.7, by = 0.03))) {
#   for (baseline_case_hosp_frac in seq(0.0001, 0.05, by=0.002)) {
#     for (time_since in c(-10, -7, -5, 0, 5, 7, 10)) {

# Parameters that go with 092425
for (lambda_1 in c(seq(0.85, 1.3, by = 0.05))) {
  for (lambda_2 in c(seq(0.85, 1.3, by = 0.05))) {
    for (baseline_case_hosp_frac in c(seq(0.00005, 0.03, by=0.002))) {
      for (time_since in c(-7, 0, 7, 14)) {
        
        # 1 MONTH
        # if(!file.exists(paste0("lambda_1", lambda_1, "-case-hosp-frac", baseline_case_hosp_frac, "-time-since", time_since, ".csv"))) {
        #   next
        # }
        # 
        # age_group_results <- read.csv(paste0("lambda_1", lambda_1, "-case-hosp-frac", baseline_case_hosp_frac, "-time-since", time_since, ".csv"))[,-1] %>%
        #   filter(type == "severe")
        
        # 2 MONTHS
        if(!file.exists(paste0("lambda_1", lambda_1, "lambda_2", lambda_2, "-case-hosp-frac", baseline_case_hosp_frac, "-time-since", time_since, ".csv"))) {
          next
        }

        age_group_results <- read.csv(paste0("lambda_1", lambda_1, "lambda_2", lambda_2, "-case-hosp-frac", baseline_case_hosp_frac, "-time-since", time_since, ".csv"))[,-1] %>%
          filter(type == "severe")
        
        combined <- rbind(age_group_results %>% 
                            filter(type == "severe") %>% 
                            dplyr::select(-type, -cases) %>%
                            mutate(group = "Simulated"),
                          observed_severe_incidence_overall %>% 
                            filter(weeks_since %in% c(0:79)) %>% 
                            rename(weeks = weeks_since, inc = observed_inc) %>% 
                            dplyr::select(age_group, weeks, inc) %>% 
                            mutate(group = "Observed"))
        
        p1 <- ggplot(data = combined, aes(x = weeks, y = inc, color = age_group, linetype = group)) + # for monthly, divide weeks by 4.34
          geom_point() + geom_line() +
          ylab("Severe Incidence (per 100,000 persons)") +
          xlab("Time (weeks)") +
          ylim(0, 30) +
          xlim(1, 8) + 
          #scale_x_continuous(breaks = c(1:18))+
          ggtitle(paste0(#"Comparison between simulated vs. observed incidence\nWaning Curve: Mean\nlamba0: ", 
            lambda_1, ", CHF:", baseline_case_hosp_frac, ", time since:", time_since)) +
          theme(title = element_text(size = 10)) 
        
        plot(p1)
      }
    }
  }
}



dev.off()
