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

#Set baseline case-hospitalization fraction and time-since shift
baseline_case_hosp_frac <- 0.005
time_since <- --7


#Load in model inputs and set age-specific lambdas
# NOTE: MAKE SURE TO SELECT HISTORICAL VS. OPTIMISTIC
# NOTE: You might need to fix file paths for the script in order for the `here` function call to work
setwd(here::here())
source(here::here("analysis/4-main analysis/simulation-inputs-historical.R"))
#source(here::here("analysis/4-main analysis/simulation-inputs-optimistic.R"))


#Calibrated lambda parameters (loaded in)
lambda_1 <- 1
lambda_2 <- 0.999239300669564
lambda_3 <- 1.00131677474849
lambda_4 <- 1.00171917687992
lambda_5 <- 1.00865700800111
lambda_6 <- 1.00029693936008
lambda_7 <- 0.997286977916913
lambda_8 <- 0.996986977916913
lambda_9 <- 1.00400330011814
lambda_10 <- 0.997251454980608
lambda_11 <- 1.001
lambda_12 <- 0.997204999183731
lambda_13 <- 1.00023367016819
lambda_14 <- 1.00331336257227
lambda_15 <- 0.997423744845859
lambda_16 <- 1.00302936138868 
lambda_17 <- 1.001
lambda_18 <- 1.0016


#Call scripts
# NOTE: You might need to fix file paths for the script in order for the `here` function call to work
source(here::here("vaccine-uptake-assignment-scenarios.R"))
source(here::here("simulation-functions-updated.R"))


#Create directory to save results
# NOTE: If you running only the main analysis under `realistic` or `optimistic` coverage, make sure to specify the file path in the write.csv line for each job (inside `if` statements below)
dir.create("final-analysis/simulation-results")

for (pop_strat in c("realistic", "optimistic")) {
  dir.create(paste0("final-analysis/simulation-results/", pop_strat))
}


if(job_selector %in% c(1:9)) {
set.seed(job_selector)
inspection_real <- simulation_semiannual_strat_8_9_10(historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
write.csv(inspection_real, paste0("final-analysis/simulation-results/strat_real-updated-", job_selector, ".csv"))
}


lambda_1 <- 1
inspection_real <- simulation_semiannual_strat_8_9_10(historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
#Average severe risk and severe protection at model start, by age and risk group
severe_risk_and_protection <- merge(protection_at_model_init[[4]], average_severe_incidence, by = c("age_group", "risk_group"), all.x = TRUE) %>%
  mutate(severe_risk = (severe_inc/100000)/7) #convert weekly severe incidence (per 100,000) into per-person daily risk

nonsevere_risk_and_protection <- merge(protection_at_model_init[[4]], average_nonsevere_incidence, by = c("age_group", "risk_group"), all.x = TRUE) %>%
  mutate(nonsevere_risk = (nonsevere_inc/100000)/7, #convert weekly nonsevere incidence (per 100,000) into per-person daily risk
         immunity_offset = (1 - mean_severe_ve)/(1 - mean_nonsevere_ve), #immunity offset by age and risk group (different from the one used for I_jt)
         nonsevere_risk_adj = nonsevere_risk / immunity_offset) #Adjust the per-person daily nonsevere risk for immunity offset at model initialization


inspection <- merge(merge(inspection_real %>% select(age_group, risk_group, total_pop, day1, nonsevere_day1),
                    nonsevere_risk_and_protection %>% select(age_group, risk_group, nonsevere_risk, immunity_offset, nonsevere_risk_adj), by = c("age_group", "risk_group"), all.x = TRUE),
                    severe_risk_and_protection %>% select(age_group, risk_group, severe_risk), by = c("age_group", "risk_group"), all.x = TRUE) %>%
  mutate(total_expected_severe = if_else(age_group == "0-17 years", 0, severe_risk * total_pop),
         total_expected_nonsevere = nonsevere_risk * total_pop,
         total_expected_nonsevere_adj = total_expected_nonsevere * 1/immunity_offset) %>%
  rename(simulated_severe = day1, simulated_nonsevere = nonsevere_day1) %>%
  select(age_group, risk_group, total_pop, total_expected_severe, total_expected_nonsevere, immunity_offset, total_expected_nonsevere_adj, simulated_severe, simulated_nonsevere)
  

lambda_1 <- 1
inspection_real <- simulation_semiannual_strat_8_9_10(historical_vax_assignment(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))


age_group_summarised <- inspection_real %>% group_by(age_group) %>% summarise(across(c(-"risk_group"), sum)) 
age_group_results <- merge(melt(age_group_summarised, id.vars = c("age_group")) %>% filter(variable != "total_pop") %>%
                             mutate(day = readr::parse_number(as.character(variable)),
                                    type = if_else(str_detect(variable, "nonsevere"), "nonsevere", "severe"),
                                    weeks = ceiling(day/7)) %>% 
                             group_by(weeks, type, age_group) %>% 
                             summarise(value = sum(value)) ,
                           melt(age_group_summarised, id.vars = c("age_group")) %>% filter(variable == "total_pop") %>% dplyr::select(-variable),
                           by = "age_group", all.x = TRUE) %>%
  mutate(inc = value.x/value.y * 100000) %>% dplyr::select(-c("value.x", "value.y")) %>%
  mutate(inc = if_else(weeks == 79, inc * 7, inc))


combined <- rbind(age_group_results %>% filter(type == "severe") %>% dplyr::select(-type) %>% mutate(group = "Simulated"),
                  weekly_severe_incidence %>% filter(weeks_since %in% c(0:79)) %>% rename(weeks = weeks_since,
                                                                                          inc = adj_inc) %>% dplyr::select(age_group, weeks, inc) %>% mutate(group = "Observed"))

ggplot(data = combined, aes(x = weeks/4.34, y = inc, color = age_group, linetype = group)) +
  geom_point() + geom_line() +
  #geom_vline(xintercept = 52) +
  ylab("Severe Incidence (per 100,000 persons)") +
  xlab("Time (months)") +
  ylim(0, NA) +
  xlim(0, 18) + 
  #ggtitle("Simulated Severe Incidence") +
  scale_x_continuous(breaks = c(1:18))+
  ggtitle("Comparison between simulated vs. observed incidence\nWaning Curve: Mean")+
  theme(title = element_text(size = 10)) 


if (job_selector %in% c(10:19)) {
  set.seed(job_selector %% 10)
  inspection_0 <- simulation_annual(strategy_0(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_0, paste0("final-analysis/simulation-results/realistic/strat_0-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(20:29)) {
  set.seed(job_selector %% 10)
  inspection_1 <- simulation_annual(strategy_1(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_1, paste0("final-analysis/simulation-results/realistic/strat_1-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(30:39)) {
  set.seed(job_selector %% 10)
  inspection_2 <- simulation_annual(strategy_2(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_2, paste0("final-analysis/simulation-results/realistic/strat_2-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(40:49)) {
  set.seed(job_selector %% 10)
  inspection_3 <- simulation_annual(strategy_3(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_3, paste0("final-analysis/simulation-results/realistic/strat_3-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(50:59)) {
  set.seed(job_selector %% 10)
  inspection_4 <- simulation_annual(strategy_4(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_4, paste0("final-analysis/simulation-results/realistic/strat_4-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(60:69)) {
  set.seed(job_selector %% 10)
  inspection_5 <- simulation_semiannual_strat_5_7(strategy_5(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_5, paste0("final-analysis/simulation-results/realistic/strat_5-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(70:79)) {
  set.seed(job_selector %% 10)
  inspection_6 <- simulation_annual(strategy_6(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_6, paste0("final-analysis/simulation-results/realistic/strat_6-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(80:89)) {
  set.seed(job_selector %% 10)
  inspection_7 <- simulation_semiannual_strat_5_7(strategy_7(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_7, paste0("final-analysis/simulation-results/realistic/strat_7-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(90:99)) {
  set.seed(job_selector %% 10)
  inspection_8 <- simulation_semiannual_strat_8_9_10(strategy_8(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_8, paste0("final-analysis/simulation-results/realistic/strat_8-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(100:109)) {
  set.seed(job_selector %% 10)
  inspection_9 <- simulation_semiannual_strat_8_9_10(strategy_9(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_9, paste0("final-analysis/simulation-results/realistic/strat_9-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(110:119)) {
  set.seed(job_selector %% 10)
  inspection_10 <- simulation_semiannual_strat_8_9_10(strategy_10(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_10, paste0("final-analysis/simulation-results/realistic/strat_10-", job_selector %% 10, ".csv"))
}


if (job_selector %in% c(120:129)) {
  set.seed(job_selector %% 10)
  inspection_1 <- simulation_annual(strategy_1(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_1, paste0("final-analysis/simulation-results/optimistic/strat_1-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(130:139)) {
  set.seed(job_selector %% 10)
  inspection_2 <- simulation_annual(strategy_2(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_2, paste0("final-analysis/simulation-results/optimistic/strat_2-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(140:149)) {
  set.seed(job_selector %% 10)
  inspection_3 <- simulation_annual(strategy_3(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_3, paste0("final-analysis/simulation-results/optimistic/strat_3-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(150:159)) {
  set.seed(job_selector %% 10)
  inspection_4 <- simulation_annual(strategy_4(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_4, paste0("final-analysis/simulation-results/optimistic/strat_4-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(160:169)) {
  set.seed(job_selector %% 10)
  inspection_5 <- simulation_semiannual_strat_5_7(strategy_5(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_5, paste0("final-analysis/simulation-results/optimistic/strat_5-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(170:179)) {
  set.seed(job_selector %% 10)
  inspection_6 <- simulation_annual(strategy_6(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_6, paste0("final-analysis/simulation-results/optimistic/strat_6-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(180:189)) {
  set.seed(job_selector %% 10)
  inspection_7 <- simulation_semiannual_strat_5_7(strategy_7(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_7, paste0("final-analysis/simulation-results/optimistic/strat_7-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(190:199)) {
  set.seed(job_selector %% 10)
  inspection_8 <- simulation_semiannual_strat_8_9_10(strategy_8(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_8, paste0("final-analysis/simulation-results/optimistic/strat_8-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(200:209)) {
  set.seed(job_selector %% 10)
  inspection_9 <- simulation_semiannual_strat_8_9_10(strategy_9(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_9, paste0("final-analysis/simulation-results/optimistic/strat_9-", job_selector %% 10, ".csv"))
}

if (job_selector %in% c(210:219)) {
  set.seed(job_selector %% 10)
  inspection_10 <- simulation_semiannual_strat_8_9_10(strategy_10(clean_df[[1]], first_dose_coverage, second_dose_coverage, next_year_dose_data, realistic_ind))
  write.csv(inspection_10, paste0("final-analysis/simulation-results/optimistic/strat_10-", job_selector %% 10, ".csv"))
}

