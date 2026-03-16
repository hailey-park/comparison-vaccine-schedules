###################################################################################################
#Title: COVID-NET severe COVID-19 incidence cleaning 
#Author: Hailey Park
#Date: December 4, 2024
###################################################################################################

rm(list=ls())

#setwd(here::here())

#Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)
library(scales)
library(data.table)
library(gridExtra)

setwd("~/comparison-vaccine-schedules/")
# Load in data
# old versions: 20250321
# more recent versions: 20260218
covidnet_char_data <- read.csv("data/raw-data/Patient_Characteristics_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20260218.csv")
covidnet_data <- read.csv("data/raw-data/Weekly_Rates_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20260218.csv")
prevalence_data <- read.csv("data/processed-data/risk_group_prevalence_estimates_1+comorbidity.csv")

#Estimating monthly incidence of COVID-19 hospitalizations by age group and risk, adjusting for COVID-19-reason hospitalizations
#NOTE:    Used data from CDC COVID-NET Surveillance system on COVID-19 Hospitalizations.
#         Link: https://data.cdc.gov/Public-Health-Surveillance/Weekly-Rates-of-Laboratory-Confirmed-COVID-19-Hosp/6jg4-xsqq/about_data
#         Link: https://data.cdc.gov/Public-Health-Surveillance/Patient-Characteristics-of-Laboratory-Confirmed-CO/bigw-pgk2/about_data 

#Pull relevant patient characteristics (reason admitted, underlying medical conditions)
covidnet_char_clean <- covidnet_char_data %>% 
  separate(Time, c("month", "year"), sep = "_") %>%
  filter(Time.Period == "Month",
         Sex == "Overall",
         Race.Ethnicity == "Overall",
         year %in% c("2023", "2024", "2025"),
         Estimate.Type == "Percent") %>%
  dplyr::select(Strata, Age.Category, COVID, month, year, Estimate) %>%
  mutate(month = case_when(month == "Novermber" ~ "November", 
                           month == "Febuary" ~ "February",
                           TRUE ~ month))

hosp_reason_data <- covidnet_char_clean %>% 
  filter(Strata == "Admitted for COVID") %>% 
  rename(hosp_reason_estimate = Estimate,
         age_group = Age.Category) %>%
  dplyr::select(age_group, month, year, hosp_reason_estimate) %>%
  mutate(hosp_reason_estimate = as.numeric(hosp_reason_estimate))

underlying_cond_data <- covidnet_char_clean %>% 
  filter(Strata == "Any underlying condition") %>% 
  rename(underlying_cond_estimate = Estimate,
         age_group = Age.Category) %>%
  dplyr::select(age_group, month, year, underlying_cond_estimate) %>%
  mutate(underlying_cond_estimate = as.numeric(underlying_cond_estimate))

combined_char_data <- merge(hosp_reason_data, underlying_cond_data, by = c("age_group", "month", "year"), all.x = TRUE)

#Split the 65+ age group into a 65-74 and 75+ age group
#below65 <- combined_char_data %>% filter(age_group %in% "≥65 years")
above65_74 <- combined_char_data %>% 
  filter(age_group == "≥65 years") %>%
  mutate(age_group = "65-74 years")

above75 <- combined_char_data %>% 
  filter(age_group == "≥65 years") %>%
  mutate(age_group = "≥75 years")

#Split the 18-49 year age group into 18-29 and 30-49
not_altered_age_groups <- combined_char_data %>% 
  filter(!age_group %in% c("18-49 years", "≥65 years"))

age_18_29 <- combined_char_data %>% 
  filter(age_group == "18-49 years") %>%
  mutate(age_group = "18-29 years")

age_30_49 <- combined_char_data %>% 
  filter(age_group == "18-49 years") %>%
  mutate(age_group = "30-49 years")

combined_char_data <- as.data.frame(rbind(not_altered_age_groups, age_18_29, age_30_49, above65_74, above75))

#Clean COVID-NET severe inc data
severe_inc_clean <- covidnet_data %>% 
  mutate(week = as.Date(str_sub(Week_ending_date, 1, 10)),
         weeks_since = as.numeric(interval(as.Date('2023-07-01'), week) %/% weeks(1)),
         month = format(week, '%B'),
         year = format(week, "%Y")) %>%
  filter(State == "COVID-NET",
         Season %in% c("2019-20","2020-21","2021-22", "2022-23", "2023-24", "2024-25"),
         AgeCategory %in% c("0-17 years (Children)", "18-49 years", "50-64 years", "65-74 years", "≥75 years"),
         Sex == "All",
         Race == "All",
         week >= as.Date("2023-07-01")) %>%
  dplyr::select(week, weeks_since, month, year, AgeCategory, WeeklyRate) %>%
  rename(age_group = AgeCategory,
         inc = WeeklyRate) 

#add extended age groups
age_18_29_inc <- severe_inc_clean %>% 
  filter(age_group == "18-49 years") %>%
  mutate(age_group = "18-29 years")

age_30_49_inc <- severe_inc_clean %>% 
  filter(age_group == "18-49 years") %>%
  mutate(age_group = "30-49 years")

severe_inc_clean <- data.frame(rbind(severe_inc_clean %>% 
                                       filter(age_group != "18-49 years"), age_18_29_inc, age_30_49_inc))

# Remove September 2025 because we don't have corresponding hosp_reason_estimate
severe_inc_clean <- subset(severe_inc_clean, !(year == 2025 & month == "September"))

#Clean prevalence data (3 risk groups)
# Here prevalence is the demographics of what % of the total population is in each age-risk group
# Stored in this google doc: https://docs.google.com/document/d/1Asdh0zHqd0qrPR8VAFE5jqTb70KRfA5z3eoslgYE-ug/edit?pli=1&tab=t.0 FIXME (update ref later)
prevalence_clean <- prevalence_data %>%
  pivot_wider(names_from = risk_group,
              values_from = prevalence) %>%
  mutate(age_group_prev = healthy + immunocompromised + `higher risk`)

#Add immunocompromised incidence data (https://www.cdc.gov/acip/downloads/slides-2024-10-23-24/03-COVID-Taylor-508.pdf)
immunocompromised_data <- data.frame(age_group = prevalence_clean$age_group,
                                     perc_immuno = c(14, 16, 16, 22, 15, 15))

#Combine severe inc data with patient characteristic info
combined_inc_with_char <- merge(merge(merge(severe_inc_clean, combined_char_data, by = c("age_group", "month", "year"), all.x = TRUE),
                                      prevalence_clean, id = "age_group", all.x = TRUE),
                                immunocompromised_data, id = "age_group", all.x = TRUE)

# Merge all the datasets
# Here, note the difference in denominator:
#   underlying_cond_estimate and perc_immuno are percents with respect to the age group
#   age_group_prev and risk status (e.g., 'higher risk') are percents with respect to the total population
CF <- 1 # FIXME: add more details for corrective factor, 0.65 for _2 version

adj_severe_inc_clean <- combined_inc_with_char %>%
  mutate(adj_inc = hosp_reason_estimate/100 * inc,
         higher_risk_inc = (underlying_cond_estimate - perc_immuno)/100 * adj_inc * age_group_prev/(`higher risk` + immunocompromised*(1-CF)),
         healthy_inc = (1 - underlying_cond_estimate/100) * adj_inc * age_group_prev/healthy,
         immunocompromised_inc = perc_immuno/100 * adj_inc * age_group_prev/(immunocompromised*CF),
         age_group = case_when(age_group == "0-17 years (Children)" ~ "0-17 years",
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) %>%
  dplyr::select(age_group, week, weeks_since, inc, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc)


#Save age-specific severe COVID-19 incidence data
write.csv(adj_severe_inc_clean %>% 
            dplyr::select(week, weeks_since, age_group, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc), 
          "data/clean-data/weekly-incidence-estimates-US-validationPeriod-updated.csv") # FIXME -- remove _2 if not using CF


# Create different age/risk-specific incidence curves that use an age-specific average of % underlying conditions
# to use for the model input. We use the average instead of time-varying.
# Here, note the difference in denominator:
#   underlying_cond_estimate and perc_immuno are percents with respect to the age group
#   age_group_prev and risk status (e.g., 'higher risk') are percents with respect to the total population
adj_severe_inc_clean_model_input <- combined_inc_with_char %>%
  group_by(age_group) %>%
  mutate(underlying_cond_estimate = mean(underlying_cond_estimate)) %>%
  ungroup() %>%
  mutate(adj_inc = hosp_reason_estimate/100 * inc,
         higher_risk_inc = (underlying_cond_estimate - perc_immuno)/100 * adj_inc * age_group_prev/(`higher risk` + immunocompromised*(1-CF)),
         healthy_inc = (1 - underlying_cond_estimate/100) * adj_inc * age_group_prev/healthy,
         immunocompromised_inc = perc_immuno/100 * adj_inc * age_group_prev/(immunocompromised*CF),
         age_group = case_when(age_group == "0-17 years (Children)" ~ "0-17 years",
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) %>%
  dplyr::select(age_group, week, weeks_since, inc, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc)


#Save age-specific severe COVID-19 incidence data
write.csv(adj_severe_inc_clean_model_input %>% dplyr::select(week, weeks_since, age_group, adj_inc, healthy_inc, immunocompromised_inc, higher_risk_inc),
          "data/clean-data/weekly-incidence-estimates-US-validationPeriod-modelInput-updated.csv")


#Plot
plot_data <- adj_severe_inc_clean_model_input %>%
  pivot_longer(
    cols = -c(age_group, week, weeks_since),  # all columns except these
    names_to = "variable",
    values_to = "value"
  )

p1 <- plot_data %>% 
  filter(variable %in% c("healthy_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  ylim(0, 30) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US\nData: COVID-NET\nRisk Group: Healthy")

p2 <- plot_data %>% 
  filter(variable %in% c("higher_risk_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  ylim(0, 75) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US\nData: COVID-NET\nRisk Group: Higher risk")

p3 <- plot_data %>% 
  filter(variable %in% c("immunocompromised_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  ylim(0, 100) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  labs(color = "Age Group", linetype = "Risk Group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US\nData: COVID-NET\nRisk Group: Immunocompromised")

grid.arrange(p1, p2, p3)

###################################################################################################
plot_data %>% 
  filter(variable %in% c("adj_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  #ylim(0, 100) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US")

plot_data %>% 
  filter(variable %in% c("adj_inc")) %>%
  ggplot(aes(week, value, color = age_group, linetype = variable)) + 
  geom_line() + geom_point() +
  #ylim(0, 100) +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Weekly Incidence of COVID Hospitalizations in US") + 
  facet_wrap(~age_group, scales = "free_y")

# Create a data frame with age groups and their populations
pop_df <- data.frame(
  age_group = c("0-17 years", "18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),  # or however your age groups are named
  population = c(0.2176, 0.1570, 0.2621, 0.1868, 0.1035, 0.073)
)

plot_data <- plot_data %>%
  left_join(pop_df, by = "age_group")

# Calculate population-weighted total incidence
total_incidence <- plot_data %>%
  filter(variable %in% c("adj_inc")) %>%
  group_by(week) %>%
  summarise(
    total_cases = sum(value * population / 100000),  # convert rates to actual cases
    total_population = sum(population),
    total_incidence = (total_cases / total_population) * 100000  # back to rate per 100k
  )

# Plot it
total_incidence %>%
  ggplot(aes(week, total_incidence)) + 
  geom_line() + geom_point() +
  xlab("Time") +
  ylab("Weekly incidence (per 100,000 persons)") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2023-07-01"), as.Date("2025-09-01")),
               breaks = "1 months") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Total Weekly Incidence of COVID Hospitalizations in US")
