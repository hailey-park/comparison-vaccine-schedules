###################################################################################################
#Title: Relative VE Comparison
#Author: Hailey Park
#Started: November 1, 2023
#Last updated: January 21, 2026, KMB
###################################################################################################

rm(list=ls())
setwd(here::here())

#Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)

#Relative VE Function
relative_ve_fn <- function(baseline, upper) {
  # Convert VE to relative risk
  # Take the ratio of upper/baseline relative risk
  # Convert back to VE
  
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  
  return(relative_ve)
}

######################################################################################################################################################################################################
#Read in waning relative VE literature estimates against severe COVID-19 infection
#     (read in specific literature reference and offset value)

#NOTE: These estimates are from Lin et al. (New England Journal of Medicine, 2023)
#      Offset: 7 months
#      DOI: 10.1056/NEJMc2215471; Figure S2 Panel B
#      Link: https://www.nejm.org/doi/full/10.1056/NEJMc2215471
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:5), 2),
                                   relative_ve = c(0.6471, 0.4747, 0.4495, 0.4242, 0.3951, 0.3698, 0.771, 0.4241, 0.3872, 0.3424, 0.3016 , 0.2529),
                                   prior_inf = rep(0:1, each = 6)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  dplyr::select(-prior_inf)

offset <- 7

#NOTE: These estimates are from Link-Gelles et al. (CDC slides, 10/23/24)
#      Offset: 4 months
#      Slide 27
#      Link: https://www.cdc.gov/acip/downloads/slides-2024-10-23-24/04-COVID-Link-Gelles-508.pdf 
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:10), 2),
                                   relative_ve = rep(c(0.5, 0.5, 0.5, 0.38, 0.38, 0.21, 0.21, 0, 0 ,0, 0), 2),
                                   prior_inf = rep(0:1, each = 11)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  dplyr::select(-prior_inf)

offset <- 4

#NOTE: These estimates are from Link-Gelles et al. (CDC MMWR, 2023)
#      Offset: 3 months
#      Table 2, Hospitalization, 18+ years
#      Link: https://www.cdc.gov/mmwr/volumes/72/wr/mm7221a3.htm
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:6), 2),
                                   relative_ve = rep(c(0.62, 0.62, 0.62, 0.47, 0.47, 0.47, 0.24), 2),
                                   prior_inf = rep(0:1, each = 7)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  dplyr::select(-prior_inf)

offset <- 3

#NOTE: These estimates are for IMMUNOCOMPROMISED from Link-Gelles et al. (CDC MMWR, 2023)
#      Offset: 2 months
#      Slide 16
#      Link: https://www.cdc.gov/acip/downloads/slides-2024-10-23-24/04-COVID-Link-Gelles-508.pdf 
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:6), 2),
                                   relative_ve = rep(c(0.36, 0.36, 0.36, 0.23, 0.23, 0.23, .01), 2),
                                   prior_inf = rep(0:1, each = 7)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates", "Literature Estimates")) %>%
  dplyr::select(-prior_inf)

offset <- 2

###################################################################################################
#Read in generated waning curves (absolute VE)
#     Note: Choose either immunocompetent or immunocompromised

#Immunocompetent
waning_data_hybrid <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% 
  summarise(ve_pred = mean(ve_pred))

waning_data_vaccine <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% 
  summarise(ve_pred = mean(ve_pred))

waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

# #Immunocompromised
# waning_data_hybrid <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")[,-1]) %>%
#   group_by(age_group, estimate, prior_inf, months) %>% 
#   summarise(ve_pred = mean(ve_pred))
# 
# waning_data_vaccine <- unique(read.csv("data/clean-data/severe_waning_predictions_weekly_vaccineImmunityOnly_immunocompromised.csv")[,-1]) %>%
#   group_by(age_group, estimate, prior_inf, months) %>% 
#   summarise(ve_pred = mean(ve_pred))
# 
# waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

###################################################################################################
# CHECK if relative VE under our waning curve predictions are similar to published estimates
#     Calculate relative VE by comparing protection at each time point to protection at a 
#     baseline time that is 'offset' months later (e.g., month 2 vs month 9 if offset = 7).
#     This simulates the incremental benefit of getting a new dose vs keeping an old dose.
#     The offset value depends on which literature source is being compared (see lines 34-84)

waning_data <- waning_data %>%
  mutate(group = if_else(prior_inf == 1, "Prior Infection; Mean", "No Prior Infection; Mean"))

adj_relative_waning <- waning_data %>% 
  filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = offset * 2, default = 24)) %>% 
  dplyr::select(group, months, months_baseline, ve_pred)

waning_data <- waning_data %>% 
  filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  rename(ve_pred_baseline = ve_pred) %>% 
  dplyr::select(group, months, ve_pred_baseline)

adj_waning_clean <- merge(adj_relative_waning, waning_data, by.x = c("months_baseline", "group"),
                          by.y = c("months", "group"), all.x = TRUE) %>% 
  mutate(months = if_else(months < 1, 0.5, round(months)),
         relative_ve = relative_ve_fn(ve_pred_baseline, ve_pred),
         group = if_else(group == "Prior Infection; Mean", "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  dplyr::select(months, relative_ve, group) 

combined_rel <- rbind(adj_waning_clean, relative_ve_data_lit)

ggplot(combined_rel, aes(months, relative_ve*100, color = factor(group))) +
  geom_line(size = 0.75) +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0, 100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0, 4, 8, 12, 16, 20, 24), expand = c(0, 0)) +
  #ggtitle("Relative VE Comparison\nReference: Link-Gelles et al. (CDC Slides, 10/23/24)\nIncremental benefit of updated 2023-2024 booster\nOffset: 4 months")
  #ggtitle("Relative VE Comparison\nReference: Link-Gelles et al. (CDC MMWR, 2023)\nIncremental benefit of bivalent booster\nOffset: 3 months")
  #ggtitle("Relative VE Comparison\nReference: Lin et al. (NEJM, 2023)\nOffset: 7 months")
  ggtitle("Relative VE Comparison\nReference: Link-Gelles et al. (CDC Slides, 10/23/24)\nIncremental benefit of updated 2023-2024 booster\nOffset: 2 months\nRisk group: Immunocompromised")


######################################################################################################################################################################################################
#Checking relative waning for nonsevere curves

#NOTE: These estimates are from de Gier et al. (Nature Comms, 2023)
#      Table S1 (average of 3 and 4-event adjusted HR estimates for hybrid immunity only)
#      Link: https://www.nature.com/articles/s41467-023-40195-z#MOESM4  
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:9), 2),
                                   relative_ve = rep(c(0.78, 0.78, 0.57, 0.57, 0.57, 0.44, 0.44, 0.44, 0.15, 0.15), 2),
                                   prior_inf = rep(2, each = 10)) %>%
  mutate(group = if_else(prior_inf == 2, "Literature Estimates", NA)) %>%
  dplyr::select(-prior_inf)

offset <- 7 #not confirmed

#NOTE: These estimates are from Link-Gelles et al. (CDC MMWR, 2024)
#      Table 2, 18+ years
#      Link: https://www.cdc.gov/mmwr/volumes/73/wr/mm7304a2.htm
relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:4)),
                                   relative_ve = rep(c(0.58, 0.58, 0.58, 0.49, 0.49), 2),
                                   prior_inf = rep(2, each = 5)) %>%
  mutate(group = if_else(prior_inf == 2, "Literature Estimates", NA)) %>%
  dplyr::select(-prior_inf)

offset <- 2

###################################################################################################

#Read in waning prediction data (absolute VE)

#Immunocompetent
waning_data_hybrid <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_hybridImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data_vaccine <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_vaccineImmunityOnly.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

#Immunocompromised
waning_data_hybrid <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_hybridImmunityOnly_immunocompromised.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data_vaccine <- unique(read.csv("data/clean-data/nonsevere_waning_predictions_weekly_vaccineImmunityOnly_immunocompromised.csv")[,-1]) %>%
  group_by(age_group, estimate, prior_inf, months) %>% summarise(ve_pred = mean(ve_pred))
waning_data <- as.data.frame(rbind(waning_data_hybrid, waning_data_vaccine))

#We are checking if the relative VE under our waning curve predictions are similar to published estimates.
#We are checking the relative VE given that the average time since last immune event (dose; infection) is around 7-8 months.

waning_data <- waning_data %>%
  mutate(group = if_else(prior_inf == 1, "Prior Infection; Mean", "No Prior Infection; Mean"))

adj_relative_waning <- waning_data %>% filter(age_group == "50-64 years", estimate == 'mean') %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = offset * 2, default = 24)) %>% 
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

combined_rel <- rbind(adj_waning_clean, relative_ve_data_lit)

ggplot(combined_rel, aes(months, relative_ve*100, color = factor(group))) +
  geom_line(size = 0.75) +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0, 4, 8, 12, 16, 20, 24), expand = c(0, 0)) +
  ggtitle("Relative VE Comparison (symptomatic infection)\nReference: Link-Gelles et al. (CDC MMWR, 2024)\nOffset: 2 months")
#ggtitle("Relative VE Comparison\nReference: de Geir et al. (Nature Comms, 2023)\nOffset: 7 months (not confirmed)")
