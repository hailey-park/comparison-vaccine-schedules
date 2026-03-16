########################################################################################################################
#Title: Figures for technical appendix
#Author: Hailey Park
#Date: September 23, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(reshape2)
library(camcorder)
library(gridExtra)

########################################################################################################################
#Figure S1-S2 (Waning curves by COVID-19 outcome, immunity type, and immunocompromised status)

  #Read in data
waning_data_mean <- read.csv("data/clean-data/waning_data_clean.csv")[,-1]
waning_data_upper <- read.csv("data/clean-data/waning_data_clean_upper.csv")[,-1]
waning_data_lower <- read.csv("data/clean-data/waning_data_clean_lower.csv")[,-1]

  #Plots with 95% UI by risk group

#combine datasets
combined <- rbind(waning_data_mean %>% mutate(estimate = 'Mean'),
                  waning_data_upper %>% mutate(estimate = 'Upper 95% UI'),
                  waning_data_lower %>% mutate(estimate = 'Lower 95% UI')) %>%
  mutate(group = if_else(immuno == 1, paste0("Immunocompromised; ", estimate),
                         paste0("Immunocompetent; ", estimate)))


#Loop through each COVID-19 outcome and immunity type and create plots
plot_list = list()
groups <- c("severe_inf_ve", "severe_vacc_ve", "severe_hybrid_ve", "nonsevere_inf_ve", "nonsevere_vacc_ve", "nonsevere_hybrid_ve")

for (i in 1:6) {
  specific_group <- groups[i]
  group_data <- combined[, c(specific_group, "weeks", "group")] %>% 
    rename(pe = specific_group)

  
  p <- ggplot(group_data, aes(x = weeks/4.345, y = pe * 100, color = group)) +
    geom_line(linewidth = .75) +
    ggtitle(paste0("Protective Effectiveness Waning\n", specific_group)) +
    ylab("Protective Effectiveness (%)")+
    xlab("Time (months)") +
    labs(color='Group') +
    ylim(0,100) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=12))+
    scale_x_continuous(limits = c(0, 24), breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
    scale_color_manual(values = c(
      "lightskyblue1",
      "dodgerblue1",
      "royalblue3",
      "tomato",
      "red",
      "red4"))
  
  plot_list[[i]] <- p
  }


pdf("waning-plots.pdf")
for (i in 1:6) {
  print(plot_list[[i]])
}
dev.off()

########################################################################################################################
#Figure S3 (Relative ve waning curves comparison, severe VE)

#Relative VE Function
relative_ve_fn <- function(baseline, upper) {
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  return(relative_ve)
}

#Read in waning relative VE literature estimates
#NOTE: These estimates are from Lin et al. (New England Journal of Medicine, 2023)
#      DOI: 10.1056/NEJMc2215471; Figure S2 Panel B
#      Link: https://www.nejm.org/doi/full/10.1056/NEJMc2215471
relative_ve_data_lit_lin <- data.frame(months = rep(c(0.5, 1:5), 2),
                                       relative_ve = c(.6471, .4747, .4495, .4242, .3951, .3698, .771, .4241, .3872,.3424,.3016 ,.2529),
                                       prior_inf = rep(0:1, each = 6)) %>%
  mutate(group = if_else(prior_inf == 1, "Lin et al. (NEJM, 2023); Prior Infection", "Lin et al. (NEJM, 2023); No Prior Infection")) %>%
  dplyr::select(-prior_inf)

#NOTE: These estimates are from Link-Gelles et al. (CDC slides, 10/23/24)
#      Slide 27
#      Link: https://www.cdc.gov/acip/downloads/slides-2024-10-23-24/04-COVID-Link-Gelles-508.pdf 
relative_ve_data_lit_link_acip <- data.frame(months = (c(0.5, 1:10)),
                                   relative_ve = (c(.5, .5, .5, .38, .38, .21, .21, 0, 0 ,0, 0))) %>%
  mutate(group = "Link-Gelles et al. (CDC MMWR/ACIP slides, 2024)") 

#NOTE: These estimates are from Link-Gelles et al. (CDC MMWR, 2023)
#      Table 2, Hospitalization, 18+ years
#      Link: https://www.cdc.gov/mmwr/volumes/72/wr/mm7221a3.htm
relative_ve_data_lit_link_mmwr <- data.frame(months = (c(0.5, 1:6)),
                                   relative_ve = (c(.62, .62, .62, .47, .47, .47, .24))) %>%
  mutate(group = "Link-Gelles et al. (CDC MMWR, 2023)") 

combined_lit_data <- rbind(relative_ve_data_lit_lin, relative_ve_data_lit_link_mmwr, relative_ve_data_lit_link_acip)

#Read in waning data
waning_data_hybrid <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data_vaccine <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

#Set offset (average time since last immune events)
offset <- 7

#We are checking if the relative VE under our waning curve predictions are similar to published estimates.
#We are checking the relative VE given that the average time since last immune event (dose; infection) is around 7-8 months.
waning_data <- waning_data %>%
  mutate(group = if_else(prior_inf == 1, "Prior Infection; Mean", "No Prior Infection; Mean"))

adj_relative_waning <- waning_data %>% filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = offset * 2, default = 23)) %>% 
  dplyr::select(group, months, months_baseline, ve_pred)

waning_data <- waning_data %>% filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  rename(ve_pred_baseline = ve_pred) %>% dplyr::select(group, months, ve_pred_baseline)

adj_waning_clean <- merge(adj_relative_waning, waning_data, by.x = c("months_baseline", "group"),
                          by.y = c("months", "group"), all.x = TRUE) %>% 
  mutate(months = if_else(months < 1, 0.5, round(months)),
         relative_ve = relative_ve_fn(ve_pred_baseline, ve_pred),
         group = if_else(group == "Prior Infection; Mean", "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  dplyr::select(months, relative_ve, group) 

combined_rel <- rbind(adj_waning_clean, combined_lit_data)

ggplot(combined_rel, aes(months, relative_ve*100, color = factor(group))) +
  geom_line(size = .75) +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(limits = c(0, 24), breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) 

########################################################################################################################
#Figure S4 (Relative ve waning curves comparison, nonsevere VE)

#Relative VE Function
relative_ve_fn <- function(baseline, upper) {
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  return(relative_ve)
}

#Read in waning relative VE literature estimates
#NOTE: These estimates are from de Gier et al. (Nature Comms, 2023)
#      Table S1 (average of 3 and 4-event adjusted HR estimates for hybrid immunity only)
#      Link: https://www.nature.com/articles/s41467-023-40195-z#MOESM4  
relative_ve_data_lit_gier <- data.frame(months = (c(0.5, 1:9)),
                                   relative_ve = (c(.78, .78, .57, .57, .57, .44, .44, .44, .15, .15)))  %>%
  mutate(group = "de Guir et al. (Nature Comms, 2023); Prior Infection") 

#NOTE: These estimates are from Link-Gelles et al. (CDC MMWR, 2024)
#      Table 2, 18+ years
#      Link: https://www.cdc.gov/mmwr/volumes/73/wr/mm7304a2.htm
relative_ve_data_lit_link_mmwr <- data.frame(months = (c(0.5, 1:4)),
                                   relative_ve = (c(.58, .58, .58, .49, .49))) %>%
  mutate(group = "Link-Gelles et al. (CDC MMWR, 2024)") 

combined_lit_data <- rbind(relative_ve_data_lit_gier, relative_ve_data_lit_link_mmwr)

#Read in waning data
waning_data_hybrid <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% 
  summarise(ve_pred = mean(ve_pred))
waning_data_vaccine <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% 
  summarise(ve_pred = mean(ve_pred))
waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

#Set offset (average time since last immune events)
offset <- 7

#We are checking if the relative VE under our waning curve predictions are similar to published estimates.
#We are checking the relative VE given that the average time since last immune event (dose; infection) is around 7-8 months.
waning_data <- waning_data %>%
  mutate(group = if_else(prior_inf == 1, "Prior Infection; Mean", "No Prior Infection; Mean"))

adj_relative_waning <- waning_data %>% filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = offset * 2, default = 23)) %>% 
  dplyr::select(group, months, months_baseline, ve_pred)

waning_data <- waning_data %>% filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  rename(ve_pred_baseline = ve_pred) %>% dplyr::select(group, months, ve_pred_baseline)

adj_waning_clean <- merge(adj_relative_waning, waning_data, by.x = c("months_baseline", "group"),
                          by.y = c("months", "group"), all.x = TRUE) %>% 
  mutate(months = if_else(months < 1, 0.5, round(months)),
         relative_ve = relative_ve_fn(ve_pred_baseline, ve_pred),
         group = if_else(group == "Prior Infection; Mean", "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  dplyr::select(months, relative_ve, group) 

combined_rel <- rbind(adj_waning_clean, combined_lit_data)

ggplot(combined_rel, aes(months, relative_ve*100, color = factor(group))) +
  geom_line(size = .75) +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(limits = c(0, 24), breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) 

########################################################################################################################
#Figure S5 (Time-since plots, by age group)

simulated_pop <- read.csv("data/clean-data/entire_population_model_initialization_daily_10mil.csv")[,-1] %>% 
  dplyr::select(age_group, prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) 

  
inspection <- simulated_pop %>% dplyr::select(age_group, prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
    group_by(age_group, prior_vacc, time_since_last_dose) %>% summarise(total = n()) %>% filter(prior_vacc != 'unvax')
  
  
  #plot simulated prior vacc
ggplot(data = inspection, aes(x = as.Date(time_since_last_dose), y = total * 5, color = age_group, alpha = prior_vacc, linetype = prior_vacc, group = interaction(prior_vacc, age_group))) +
  geom_line(size = 0.7) +
  ylab("Total (individuals)") +
  xlab("Time") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-12-01"), as.Date("2023-06-30")),
               breaks = "2 months") +
  scale_alpha_manual(values = c(1, 0.7, 0.7), labels = c("Boosted (1st dose)", "Boosted (2+ doses)", "Primary Series Completed"), guide = "none") +
  scale_linetype_manual(values = c(1, 2, 9), labels = c("Boosted (1st dose)", "Boosted (2+ doses)", "Primary Series Completed")) +
  scale_color_brewer(palette = "Set1") +
  ylim(0, NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(color = "Age Group",
       linetype = "Vaccination Status") +
  ggtitle(paste0("Simulated Time-Since Last Vaccination"))
  

inspection <- simulated_pop %>% dplyr::select(age_group, prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
  group_by(age_group, time_since_last_inf) %>% summarise(total = n()) %>% filter(!is.na(time_since_last_inf))


  ## plot simulated prior inf
ggplot(data = inspection, aes(x = as.Date(time_since_last_inf), y = total * 5, color = age_group, group = age_group))  +
         geom_line(size = 0.7) +
         ylab("Total (individuals)") +
         xlab("Time") +
  scale_color_brewer(palette = "Set1") +
         scale_x_date(date_labels = ("%b-%Y"),
                      limits = c(as.Date("2020-12-01"), as.Date("2023-07-01")),
                      breaks = "2 months") +
         ylim(0, 150000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(color = "Age Group") +
  ggtitle(paste0("Simulated Time-Since Last Infection"))
  
inspection <- simulated_pop %>% dplyr::select(age_group, prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
  group_by(age_group, time_since_last_reinf) %>% summarise(total = n()) %>% filter(!is.na(time_since_last_reinf))

  
  ## plot simulated prior reinf
ggplot(data = inspection, aes(x = as.Date(time_since_last_reinf), y = total * 5, color = age_group, group = age_group))  +
  geom_line(size = 0.7) +
  ylab("Total (individuals)") +
  xlab("Time") +
  scale_color_brewer(palette = "Set1") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-12-01"), as.Date("2023-07-01")),
               breaks = "2 months") +
  ylim(0, 150000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(color = "Age Group") +
  ggtitle(paste0("Simulated Time-Since Last Infection"))

 
inspection <- simulated_pop %>% dplyr::select(age_group, prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_reinf, time_since_last_dose_inf) %>%
  group_by(age_group, time_since_last_dose_inf) %>% summarise(total = n()) %>%
  mutate(time_since_adjusted = floor_date(as.Date(time_since_last_dose_inf), unit = "day")) %>%
  group_by(age_group, time_since_adjusted) %>% summarise(total = sum(total))


  
  #plot simulated time since
  
  ggplot(data = inspection, aes(x = as.Date(time_since_adjusted), y = total * 5, color = age_group, group = age_group))  +
    geom_line(size = 0.7) +
    ylab("Total (individuals)") +
    xlab("Time") +
    scale_color_brewer(palette = "Set1") +
    scale_x_date(date_labels = ("%b-%Y"),
                 limits = c(as.Date("2020-12-01"), as.Date("2023-07-01")),
                 breaks = "2 months") +
    ylim(0, NA)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(color = "Age Group") +
    ggtitle(paste0("Simulated Time-Since Last Immune Event"))
  
########################################################################################################################
#Figure S6 (Calibration plot)

  results <- results#"[INSERT DF AFTER RUNNING SIMULATION]"
  
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
    ggtitle(paste0("Comparison between simulated vs. observed incidence")) +
    theme(title = element_text(size = 10)) 
  
########################################################################################################################
#Figure S7 (Validation plot)
  
  #results <- results#"[INSERT DF AFTER RUNNING SIMULATION]"
  setwd("/Users/katebubar/comparison-vaccine-schedules/")
  # NOTE: ...validationPeriod_2 contains the corrective factor. remove _2 if using original data
  weekly_severe_incidence <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod.csv")[,-1] #This is the observed COVID-NET data where we are using time-varying (month-specific) %COVID-NET hospitalizations with underlying medical conditions
  
  #setwd("~/Desktop/final-analysis-with-CF65/")
  # Read in all 10 results files
  all_results <- list()
  for (i in 1:1) {
    all_results[[i+1]] <- read.csv(paste0("simulation-results/strat_real-updated-", i, ".csv")) %>%
      mutate(sim_run = i+1)  # Add identifier for which simulation run
  }
  
  # Combine all results
  results_combined <- bind_rows(all_results)
  
  # Process all simulations together
  age_group_summarised <- results_combined %>% 
    group_by(age_group, risk_group, sim_run) %>% 
    summarise(across(everything(), sum), .groups = "drop") 
  
  age_group_results <- merge(reshape2::melt(age_group_summarised, id.vars = c("age_group", "risk_group", "sim_run")) %>% 
                               filter(variable != "total_pop") %>%
                               filter(variable != "X") %>%
                               filter(variable != "total_vaccines")  %>%
                               filter(variable != "total_vaccinated")  %>%
                               filter(variable != "total_severe_vax")  %>%
                               filter(variable != "total_nonsevere_vax")  %>%
                               mutate(day = readr::parse_number(as.character(variable)),
                                      type = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
                                      weeks = ceiling(day/7)) %>% 
                               group_by(weeks, type, age_group, risk_group, sim_run) %>% 
                               summarise(value = sum(value), .groups = "drop"),
                             reshape2::melt(age_group_summarised, id.vars = c("age_group", "risk_group", "sim_run")) %>% 
                               filter(variable == "total_pop") %>% 
                               dplyr::select(-variable),
                             by = c("age_group", "risk_group", "sim_run"), all.x = TRUE) %>%
    mutate(inc = value.x/value.y * 100000) %>% 
    dplyr::select(-c("value.x", "value.y")) %>%
    mutate(inc = if_else(weeks == 79, inc * 7, inc))
  
  #Select which risk group to plot validation performance
  # In `age_group_results` which is the simulated data, filter `risk_group` for either "immunocompromised", "higher risk", "healthy"
  # In weekly_severe_incidence which is the observed data, set `inc` to either "immunocompromised_inc", "higher_risk_inc", or "healthy_inc"
  # In ggplot, change title to correct risk group
  # Select which risk group to plot validation performance
  
  # HEALTHY
  combined <- rbind(age_group_results %>% 
                      filter(type == "severe", risk_group == "healthy") %>% 
                      dplyr::select(-type, -risk_group) %>%
                      mutate(group = "Simulated"),
                    weekly_severe_incidence %>% 
                      filter(weeks_since %in% c(0:79)) %>% 
                      rename(weeks = weeks_since,
                             inc = healthy_inc) %>% 
                      dplyr::select(age_group, weeks, inc) %>% 
                      mutate(group = "Observed",
                             sim_run = NA))  # Add sim_run column for observed data
  
  ggplot(data = combined, aes(x = weeks/4.34, y = inc, color = age_group)) +
    geom_point(data = combined %>% filter(group == "Observed"), size = 2) + 
    geom_line(data = combined %>% filter(group == "Observed"), linewidth = 0.8) +
    geom_line(data = combined %>% filter(group == "Simulated"), 
              aes(group = interaction(age_group, sim_run)), alpha = 0.4, linewidth = 0.3) +
    ylab("Severe Incidence (per 100,000 persons)") +
    xlab("Time (months)") +
    scale_x_continuous(breaks = c(1:18))+
    ggtitle("Comparison between simulated vs. observed incidence\nRisk Group: Healthy")+
    theme(title = element_text(size = 10))
  
  # CHECK CUMULATIVE
  observed_cumul <- combined %>%
    filter(group == "Observed") %>%
    group_by(age_group) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul_by_run <- combined %>%
    filter(group == "Simulated") %>%
    group_by(age_group, sim_run) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul <- simulated_cumul_by_run %>%
    group_by(age_group) %>%
    summarize(mean_cumulative_inc = mean(cumulative_inc, na.rm = TRUE))
  
  error_data <- observed_cumul %>%
    left_join(simulated_cumul, by = "age_group") %>%
    filter(age_group != "0-17 years") %>%
    mutate(
      abs_error = abs(cumulative_inc - mean_cumulative_inc),
      rel_error = abs_error / cumulative_inc * 100,
      label = ifelse(age_group == "18-29 years",
                     paste0("AE: ", round(abs_error, 1), "\nRE: ", round(rel_error, 1), "%"),
                     paste0(round(abs_error, 1), "\n", round(rel_error, 1), "%"))
    )
  
  p1 <- ggplot() + 
    geom_point(data = simulated_cumul_by_run %>% filter(age_group != "0-17 years"), aes(x = age_group, y = cumulative_inc), color = "grey12", alpha = 0.2, size = 3) + 
    geom_point(data = simulated_cumul %>% filter(age_group != "0-17 years"), aes(x = age_group, y = mean_cumulative_inc), color = "grey12", shape = 95, size = 7) + 
    geom_point(data = observed_cumul %>% filter(age_group != "0-17 years"), aes(x = age_group, y = cumulative_inc), color = "purple", size = 2.5) +
    geom_text(data = error_data, 
              aes(x = age_group, y = max(cumulative_inc, mean_cumulative_inc) * 1.1, label = label),
              size = 3, vjust = 1) +
    theme_bw() + 
    xlab("age group")+
    ylab("cumulative hospitalizations (per 100,000)") +
    ggtitle("Healthy")
  
  # HIGHER RISK 
  combined <- rbind(age_group_results %>% 
                      filter(type == "severe", risk_group == "higher risk") %>% 
                      dplyr::select(-c(type, risk_group)) %>% 
                      mutate(group = "Simulated"),
                    weekly_severe_incidence %>% 
                      filter(weeks_since %in% c(0:79)) %>% 
                      rename(weeks = weeks_since,
                             inc = higher_risk_inc) %>% 
                      dplyr::select(age_group, weeks, inc) %>% 
                      mutate(group = "Observed",
                             sim_run = NA))
  
  ggplot(data = combined, aes(x = weeks/4.34, y = inc, color = age_group)) +
    geom_point(data = combined %>% filter(group == "Observed"), size = 2) + 
    geom_line(data = combined %>% filter(group == "Observed"), linewidth = 0.8) +
    geom_line(data = combined %>% filter(group == "Simulated"), 
              aes(group = interaction(age_group, sim_run)), alpha = 0.4, linewidth = 0.3) +
    ylab("Severe Incidence (per 100,000 persons)") +
    xlab("Time (months)") +
    xlim(0, 18) + 
    scale_x_continuous(breaks = c(1:18))+
    ggtitle("Comparison between simulated vs. observed incidence\nRisk Group: Higher risk")+
    theme(title = element_text(size = 10))
  
  # CHECK CUMULATIVE
  observed_cumul <- combined %>%
    filter(group == "Observed") %>%
    group_by(age_group) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul_by_run <- combined %>%
    filter(group == "Simulated") %>%
    group_by(age_group, sim_run) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul <- simulated_cumul_by_run %>%
    group_by(age_group) %>%
    summarize(mean_cumulative_inc = mean(cumulative_inc, na.rm = TRUE))
  
  error_data <- observed_cumul %>%
    left_join(simulated_cumul, by = "age_group") %>%
    filter(age_group != "0-17 years") %>%
    mutate(
      abs_error = abs(cumulative_inc - mean_cumulative_inc),
      rel_error = abs_error / cumulative_inc * 100,
      label = ifelse(age_group == "18-29 years",
                     paste0("AE: ", round(abs_error, 1), "\nRE: ", round(rel_error, 1), "%"),
                     paste0(round(abs_error, 1), "\n", round(rel_error, 1), "%"))
    )
  
  p2 <- ggplot() + 
    geom_point(data = simulated_cumul_by_run %>% filter(age_group != "0-17 years"), aes(x = age_group, y = cumulative_inc), color = "grey12", alpha = 0.2, size = 3) + 
    geom_point(data = simulated_cumul %>% filter(age_group != "0-17 years"), aes(x = age_group, y = mean_cumulative_inc), color = "grey12", shape = 95, size = 7) + 
    geom_point(data = observed_cumul %>% filter(age_group != "0-17 years"), aes(x = age_group, y = cumulative_inc), color = "purple", size = 2.5) +
    theme_bw() + 
    geom_text(data = error_data, 
              aes(x = age_group, y = max(cumulative_inc, mean_cumulative_inc) * 1.1, label = label),
              size = 3, vjust = 1) +
    xlab("age group")+
    ylab("cumulative hospitalizations (per 100,000)") +
    ggtitle("Higher risk")
  
  # IMMUNOCOMPROMISED
  combined <- rbind(age_group_results %>% 
                      filter(type == "severe", risk_group == "immunocompromised") %>% 
                      dplyr::select(-c(type, risk_group)) %>% 
                      mutate(group = "Simulated"),
                    weekly_severe_incidence %>% 
                      filter(weeks_since %in% c(0:79)) %>% 
                      rename(weeks = weeks_since,
                             inc = immunocompromised_inc) %>% 
                      dplyr::select(age_group, weeks, inc) %>% 
                      mutate(group = "Observed", 
                             sim_run = NA))
  
  ggplot(data = combined, aes(x = weeks/4.34, y = inc, color = age_group)) +
    geom_point(data = combined %>% filter(group == "Observed"), size = 2) + 
    geom_line(data = combined %>% filter(group == "Observed"), linewidth = 0.8) +
    geom_line(data = combined %>% filter(group == "Simulated"), 
              aes(group = interaction(age_group, sim_run)), alpha = 0.7, linewidth = 0.3) +
    ylab("Severe Incidence (per 100,000 persons)") +
    xlab("Time (months)") +
    xlim(0, 18) + 
    scale_x_continuous(breaks = c(1:18))+
    ggtitle("Comparison between simulated vs. observed incidence\nRisk Group: Immunocompromised") +
    theme(title = element_text(size = 10))

  # CHECK CUMULATIVE
  observed_cumul <- combined %>%
    filter(group == "Observed") %>%
    group_by(age_group) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul_by_run <- combined %>%
    filter(group == "Simulated") %>%
    group_by(age_group, sim_run) %>%
    summarize(cumulative_inc = sum(inc, na.rm = TRUE))
  
  simulated_cumul <- simulated_cumul_by_run %>%
    group_by(age_group) %>%
    summarize(mean_cumulative_inc = mean(cumulative_inc, na.rm = TRUE))
  
  error_data <- observed_cumul %>%
    left_join(simulated_cumul, by = "age_group") %>%
    filter(age_group != "0-17 years") %>%
    mutate(
      abs_error = abs(cumulative_inc - mean_cumulative_inc),
      rel_error = abs_error / cumulative_inc * 100,
      label = ifelse(age_group == "18-29 years",
                     paste0("AE: ", round(abs_error, 1), "\nRE: ", round(rel_error, 1), "%"),
                     paste0(round(abs_error, 1), "\n", round(rel_error, 1), "%"))
    )
  
  p3 <- ggplot() + 
    geom_point(data = simulated_cumul_by_run %>% filter(age_group != "0-17 years"), 
               aes(x = age_group, y = cumulative_inc), color = "grey12", alpha = 0.2, size = 3) + 
    geom_point(data = simulated_cumul %>% filter(age_group != "0-17 years"), 
               aes(x = age_group, y = mean_cumulative_inc), color = "grey12", shape = 95, size = 7) + 
    geom_point(data = observed_cumul %>% filter(age_group != "0-17 years"), 
               aes(x = age_group, y = cumulative_inc), color = "purple", size = 2.5) +
    geom_text(data = error_data, 
              aes(x = age_group, y = max(cumulative_inc, mean_cumulative_inc) * 1.1, label = label),
              size = 3, vjust = 1) +  # adjust vjust to move text up/down
    theme_bw() + 
    xlab("age group")+
    ylab("cumulative hospitalizations (per 100,000)") +
    ggtitle("Immunocompromised")
  grid.arrange(p1, p2, p3, ncol = 3)
  
  
########################################################################################################################
#Figure S8-S9 (Sensitivity analysis waning curves...TBD)

