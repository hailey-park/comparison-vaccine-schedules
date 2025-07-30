########################################################################################################################
#Title: Figures for primary outcomes (vaccination strategy comparison)
#Author: Hailey Park
#Date: April 20, 2025
########################################################################################################################

rm(list=ls())

#Loading in libraries
library(tidyverse)
library(reshape2)
library(data.table)

#Read in strategies
strat_real <- read.csv("simulation-results-7/strat_real-updated.csv")[,-1]
strat_0 <- read.csv("simulation-results-7/pessimistic/strat_0.csv")[,-1]
strat_1 <- read.csv("simulation-results-7/pessimistic/strat_1.csv")[,-1]
strat_2 <- read.csv("simulation-results-7/pessimistic/strat_2.csv")[,-1]
strat_3 <- read.csv("simulation-results-7/pessimistic/strat_3.csv")[,-1]
strat_4 <- read.csv("simulation-results-7/pessimistic/strat_4.csv")[,-1]
strat_5 <- read.csv("simulation-results-7/pessimistic/strat_5.csv")[,-1]
strat_6 <- read.csv("simulation-results-7/pessimistic/strat_6.csv")[,-1]
strat_7 <- read.csv("simulation-results-7/pessimistic/strat_7.csv")[,-1]
strat_8 <- read.csv("simulation-results-7/pessimistic/strat_8.csv")[,-1]
strat_9 <- read.csv("simulation-results-7/pessimistic/strat_9.csv")[,-1]
strat_10 <- read.csv("simulation-results-7/pessimistic/strat_10.csv")[,-1]

#cleaning results (general risk group-specific)
cleaning_sim_general <- function(df) {
  clean_df <- df %>% mutate(severe_total = select(., day1:day547) %>% rowSums(na.rm = TRUE),
                            nonsevere_total = select(., nonsevere_day1:nonsevere_day547) %>% rowSums(na.rm = TRUE)) %>%
    dplyr::select(age_group, risk_group, total_pop, severe_total, nonsevere_total, total_vaccines, total_vaccinated, total_severe_vax, total_nonsevere_vax) %>%
    mutate(general_risk_cat = case_when(age_group %in% c("65-74 years", "75+ years") ~ "65+ years",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "healthy" ~ "18-64 years, healthy",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "higher risk" ~ "18-64 years, higher risk",
                                        age_group %in% c("18-29 years", "30-49 years", "50-64 years") & risk_group == "immunocompromised" ~ "18-64 years, immunocompromised",
                                        TRUE ~ age_group)) %>%
    group_by(general_risk_cat) %>% summarise(across(total_pop:total_nonsevere_vax, sum)) %>%
    mutate(severe_risk = (severe_total/total_pop)/1.5 * 100000,
           nonsevere_risk = (nonsevere_total/total_pop)/1.5 * 100000,
           severe_risk_vax = (total_severe_vax/total_vaccinated)/1.5 * 100000) %>%
    filter(general_risk_cat != "0-17 years") %>%
    dplyr::select(general_risk_cat, severe_risk, severe_risk_vax)
  
  return(clean_df)
  
}


cleaned_results <- list(strat_real, strat_0, strat_1, strat_2, strat_3, strat_4, strat_5, strat_6, strat_7, strat_8, strat_9, strat_10) %>%
  lapply(cleaning_sim_general)


combined_strategies <- rbind(cleaned_results[[1]] %>% mutate(strategy = "Historical vaccination\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[2]] %>% mutate(strategy = "No vaccination (counterfactual)\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[4]] %>% mutate(strategy = "Annual vaccine 18+ years\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[6]] %>% mutate(strategy = "Annual vaccine 65+ years,\nimmunocompromised\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[7]] %>% mutate(strategy = "Semiannual vaccine 65+ years,\nimmunocompromised\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[8]] %>% mutate(strategy = "Annual vaccine 65+ years, \nimmunocompromised + higher risk\n", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax),
                             cleaned_results[[11]] %>% mutate(strategy = "Annual vaccine 18+ years, \n2nd dose in 65+ years, immunocompromised", severe_vax_ARR = cleaned_results[[2]]$severe_risk - severe_risk_vax))


plot_data <- combined_strategies %>%
                mutate(strategy = factor(strategy, levels = c("Historical vaccination\n", 
                                                              "No vaccination (counterfactual)\n",
                                                              "Annual vaccine 18+ years\n",
                                                              "Annual vaccine 65+ years,\nimmunocompromised\n",
                                                              "Semiannual vaccine 65+ years,\nimmunocompromised\n",
                                                              "Annual vaccine 65+ years, \nimmunocompromised + higher risk\n",
                                                              "Annual vaccine 18+ years, \n2nd dose in 65+ years, immunocompromised")),
                       general_risk_cat = factor(general_risk_cat, levels = c("18-64 years, healthy",
                                                                              "18-64 years, higher risk",
                                                                              "18-64 years, immunocompromised",
                                                                              "65+ years")))
                        

ggplot(plot_data, aes(x = strategy, y = (severe_risk), color = strategy)) + 
  #facet_wrap(. ~ general_risk_cat, scales = "free", ncol = 4, strip.position = "bottom") + 
  facet_grid(. ~ general_risk_cat, switch = "both") + 
  geom_point(size = 2) +
  ylim(0, 1000) +
  scale_color_brewer(palette = "Dark2") +
  ylab("Annual Risk of Severe COVID-19 \n (cases per 100,000)") +
  xlab("General Risk Category")+
  labs(color = "Strategies")+
  #ggtitle("Annual Risk Comparison\nPopulation: everyone\nUptake: realistic coverage") +
  #ggtitle("Absolute scale, y-axis fixed") +
  theme(
    text = element_text(size=10),
    axis.text.x=element_blank())

ggplot(plot_data %>% filter(strategy != "No vaccination (counterfactual)\n"), aes(x = strategy, y = (severe_vax_ARR), color = strategy)) + 
  #facet_wrap(. ~ general_risk_cat, scales = "free", ncol = 4, strip.position = "bottom") + 
  facet_grid(. ~ general_risk_cat, switch = "both") + 
  geom_point(size = 2) +
  ylim(0, 800) +
  scale_color_manual(values = c("#1B9E77FF", "#7570B3FF", "#E7298AFF", "#66A61EFF", "#E6AB02FF", "#A6761DFF", "#666666FF")) + 
  #scale_color_brewer(palette = "Dark2") +
  ylab("Absolute Risk Averted in Vaccinated Population\n (cases per 100,000)") +
  xlab("General Risk Category")+
  labs(color = "Strategies")+
  #ggtitle("Annual Risk Comparison\nPopulation: everyone\nUptake: realistic coverage") +
  #ggtitle("Absolute scale, y-axis fixed") +
  theme(
    text = element_text(size=10),
    axis.text.x=element_blank())

