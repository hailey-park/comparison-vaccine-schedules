########################################################################################################################
#Title: Plotting script
#Author: Hailey Park
#Date: September 11, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(reshape2)
library(data.table)

set.seed(88)
results <- simulation_semiannual_strat_8_9_10(historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))

#Model init table
cleaned_results <- merge(merge(results %>% select(age_group, risk_group, total_pop, day1, nonsevere_day1),
                               average_severe_incidence, by = c("age_group", "risk_group"), all.x = TRUE),
                         average_nonsevere_incidence, by = c("age_group", "risk_group"), all.x = TRUE) %>%
  mutate(total_expected_severe = total_pop * severe_inc/100000 / 7,
         total_expected_nonsevere = total_pop * nonsevere_inc/100000 / 7) %>%
  rename(simulated_severe_inf = day1,
         simulated_nonsevere_inf = nonsevere_day1) %>%
  select(age_group, risk_group, total_pop, severe_inc, nonsevere_inc, total_expected_severe, total_expected_nonsevere, simulated_severe_inf, simulated_nonsevere_inf)

#Model table averaged over risk groups
age_group_summary <- cleaned_results %>%
  group_by(age_group) %>%
  summarise(
    simulated_severe_inf = sum(simulated_severe_inf), #sum(simulated_severe_inf * total_pop) / sum(total_pop),
    total_expected_severe = sum(total_expected_severe), #sum(total_expected_severe * total_pop) / sum(total_pop),
    total_pop = sum(total_pop),
    .groups = "drop"
  )


#results <- results#"[INSERT DF AFTER RUNNING SIMULATION]"

age_group_summarised <- results %>% 
  group_by(age_group) %>% 
  summarise(across(c(-"risk_group"), sum)) 

age_group_results <- merge(
  # incidence part
  age_group_summarised %>%
    pivot_longer(
      cols = -age_group,
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(variable != "total_pop") %>%
    mutate(
      day = parse_number(variable),
      type = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
      weeks = ceiling(day / 7)
    ) %>%
    group_by(weeks, type, age_group) %>%
    summarise(value = sum(value), .groups = "drop"),
  
  # total_pop part
  age_group_summarised %>%
    pivot_longer(
      cols = -age_group,
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(variable == "total_pop") %>%
    select(-variable),
  
  by = "age_group",
  all.x = TRUE
) %>%
  mutate(
    inc = value.x / value.y * 100000,
    inc = if_else(weeks == 79, inc * 7, inc)
  ) %>%
  select(-c("value.x", "value.y"))

combined <- rbind(age_group_results %>% 
                    filter(type == "severe") %>% 
                    dplyr::select(-type) %>%
                    mutate(group = "Simulated"),
                  weekly_severe_incidence %>% 
                    filter(weeks_since %in% c(0:79)) %>% 
                    rename(weeks = weeks_since, inc = adj_inc) %>% 
                    dplyr::select(age_group, weeks, inc) %>% 
                    mutate(group = "Observed"))

ggplot(data = combined, aes(x = weeks/4.34, y = inc, color = age_group, linetype = group)) +
  geom_point() + geom_line() +
  ylab("Severe Incidence (per 100,000 persons)") +
  xlab("Time (months)") +
  ylim(0, 60) +
  xlim(0, 18) + 
  scale_x_continuous(breaks = c(1:18))+
  ggtitle(paste0("Comparison between simulated vs. observed incidence\nWaning Curve: Mean\nlamba0: ", lambda_1, ", CHF:", baseline_case_hosp_frac, ", time since:", time_since)) +
  theme(title = element_text(size = 10)) 
