###################################################################################################
#Title: Calibration data inputs
#Author: Hailey Park
#Date: December 20, 2024
###################################################################################################

#Severe incidence estimates by age
weekly_severe_incidence <- read.csv("data/clean-data/weekly-incidence-estimates-US-validationPeriod.csv")[,-1] %>%
  mutate(age_group = case_when(age_group == "0-17 years (Children)" ~  "0-17 years", 
                               age_group == "≥75 years" ~ "75+ years",
                               TRUE ~ age_group)) 

#Ratios between age-specific case hospitalization fractions
#Baseline: 50-64 years
age_ratio_case_hosp_frac <- data.frame(age_group = c("0-17 years", "18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),
                                       case_hosp_frac = c(0.155, 0.4283, 1.545, 5.9, 13, 30)) %>%
  mutate(ratio = case_hosp_frac/5.9)

#Nonsevere incidence estimates by age (assuming same across risk groups)
#NOTE: `baseline_case_hosp_frac` will be a calibrated parameter
weekly_nonsevere_incidence <- merge(weekly_severe_incidence, age_ratio_case_hosp_frac, by = "age_group", all.x = TRUE) %>%
  mutate(healthy_inc = adj_inc * (1/(ratio * baseline_case_hosp_frac)),
         immunocompromised_inc = adj_inc * (1/(ratio * baseline_case_hosp_frac)),
         higher_risk_inc = adj_inc * (1/(ratio * baseline_case_hosp_frac))) 


#Contact matrix (POLYMOD)
contact_matrix <- read.csv("data/processed-data/contact matrix updated.csv") 

#Waning curves
waning_data_clean <- setDT(read.csv("data/clean-data/waning_data_clean.csv")[,-1])
waning_data_clean_alt <- setDT(read.csv("data/clean-data/waning_data_clean_alt.csv")[,-1]) %>%
  filter(weeks %% 2 == 0) %>% mutate(weeks = weeks/2)
setkeyv(waning_data_clean_alt, c("immuno", "weeks", "prior_inf", "prior_vacc"))

#MAKE SURE YOU ARE READING IN THE CORRECT CALIBRATION FILE
entire_pop <-  read.csv(paste0("data/clean-data/entire_population_model_initialization_daily_50mil_updated.csv"))[,-1] 

#Calculating total current infections at model initialization, based on average severe and nonsevere incidence estimates at July 2023
average_severe_incidence <- melt(weekly_severe_incidence %>% filter(weeks_since < 4) %>% group_by(age_group) %>% 
  summarise(across(healthy_inc:higher_risk_inc, mean)), id = "age_group") %>% 
  rowwise() %>% dplyr::mutate(risk_group = strsplit(as.character(variable), split="_")[[1]][1],
                              risk_group = if_else(risk_group == "higher", "higher risk", risk_group)) %>%
  rename(severe_inc = value) %>% dplyr::select(-variable)

average_nonsevere_incidence <- melt(weekly_nonsevere_incidence %>% filter(weeks_since < 4) %>% group_by(age_group) %>% 
  summarise(across(healthy_inc:higher_risk_inc, mean)), id = "age_group") %>% 
  rowwise() %>% dplyr::mutate(risk_group = strsplit(as.character(variable), split="_")[[1]][1],
                              risk_group = if_else(risk_group == "higher", "higher risk", risk_group)) %>%
  rename(nonsevere_inc = value) %>% dplyr::select(-variable)

#Estimate the total currently circulating infections @ t = 0 using the average severe and nonsevere incidences
inf_by_age <- merge(merge(entire_pop %>% group_by(age_group, risk_group) %>% summarise(total_pop = n()),
                          average_severe_incidence, by = c("age_group", "risk_group"), all.x = TRUE),
                    average_nonsevere_incidence, by = c("age_group", "risk_group"), all.x = TRUE) %>% 
  mutate(nonsevere_inf = (total_pop * nonsevere_inc/100000)/7 , #dividing by 7 to convert weekly to daily
         severe_inf = if_else(age_group == "0-17 years", 0, (total_pop * severe_inc/100000)/7 ), 
         total_inf = (nonsevere_inf + severe_inf) * 5) %>% #multiplying by 5 because infectious period is 5 days
  group_by(age_group) %>% summarise(total_inf = ceiling(sum(total_inf)),
                                    total_pop = sum(total_pop))

#Estimate factors to account for differences in magnitude after including the dynamic term
contact_matrix_adj <- data.frame(age_group = c("0-17 years","18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),
                                 contact_matrix_adj = (sum(inf_by_age$total_inf)/50000000)/c(sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X0.17.years),
                                                                                                 sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X18.29.years),
                                                                                                 sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X30.49.years),
                                                                                                 sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X50.64.years),
                                                                                                 sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X65.74.years),
                                                                                                 sum((inf_by_age$total_inf/inf_by_age$total_pop) * contact_matrix$X75..years)))

#Tracking the last 8 days of infection counts because 3 days (latent period) and 5 days (infectious period)
infection_tracker_df <- data.frame(days_since = c(1:8),
                                   age_0_17 = rep(inf_by_age$total_inf[1]/5, 8),
                                   age_18_29 = rep(inf_by_age$total_inf[2]/5, 8),
                                   age_30_49 = rep(inf_by_age$total_inf[3]/5, 8),
                                   age_50_64 = rep(inf_by_age$total_inf[4]/5, 8),
                                   age_65_74 = rep(inf_by_age$total_inf[5]/5, 8),
                                   age_75_plus = rep(inf_by_age$total_inf[6]/5, 8))
############################################################################################
#clean age matrices
clean_age_matrix <- function(df){
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
  lapply(clean_age_matrix) 

rm(entire_pop)
############################################################################################
#calculate multiplier adjustments (this accounts for the difference in magnitude between severe and nonsevere protection, by risk group)
protection_at_model_initialization <- function(df){
  with_protection <- merge((df %>% rowwise() %>%
                             mutate(weeks_since_last_dose_inf = min(max(1, weeks_since_last_dose_inf), 730))), waning_data_clean, by.x = c("immuno","weeks_since_last_dose_inf"), 
                           by.y = c("immuno", "weeks"), all.x = TRUE) %>% arrange(individual)
  
  #Individuals who are unvaccinated and no prior infection history has no protection
  immune_naive_index <- which(with_protection$prior_vacc =="unvax" & with_protection$prior_inf == 'noinf')
  with_protection[immune_naive_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- 0
  
  #Individuals with infection <3 months from simulation start has perfect immunity
  perfect_immunity_index <- which(with_protection$weeks_since_last_dose_inf < with_protection$weeks_since_last_dose & with_protection$prior_inf != "noinf"  & with_protection$weeks_since_last_dose_inf %in% c(1:13))
  with_protection[perfect_immunity_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- 1
  
  #Individuals who are unvaccinated and have prior infection history have prior infection only waning immunity
  prior_inf_only_index <- which(with_protection$prior_vacc=="unvax" & with_protection$prior_inf != "noinf")
  with_protection[prior_inf_only_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- with_protection[prior_inf_only_index, c("severe_inf_ve", "nonsevere_inf_ve")]
  
  #Individuals who are previously infected and have prior vaccination history have vaccine only waning immunity
  vaccine_only_index <- which(with_protection$prior_vacc !="unvax" & with_protection$prior_inf == "noinf")
  with_protection[vaccine_only_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- with_protection[vaccine_only_index, c("severe_vacc_ve", "nonsevere_vacc_ve")]
  
  #Individuals who have prior vaccination and have prior infection history have hybrid waning immunity
  hybrid_index <- which(with_protection$prior_vacc !="unvax" & with_protection$prior_inf != "noinf")
  with_protection[hybrid_index, c("severe_ve_pred", "nonsevere_ve_pred")] <- with_protection[hybrid_index, c("severe_hybrid_ve", "nonsevere_hybrid_ve")]
  
  group_multiplier_adj <- with_protection %>% group_by(age_group, risk_group) %>% summarise(mean_severe_ve = mean(severe_ve_pred),
                                                                                                   mean_nonsevere_ve = mean(nonsevere_ve_pred)) %>%
    mutate(multiplier_adj = (1-mean_severe_ve)/(1-mean_nonsevere_ve))
  
  
  age_group_multiplier_adj <- with_protection %>% group_by(age_group) %>% summarise(mean_severe_ve = mean(severe_ve_pred),
                                                                                    mean_nonsevere_ve = mean(nonsevere_ve_pred)) 
  
  return(list(group_multiplier_adj %>% dplyr::select(age_group, risk_group, mean_nonsevere_ve, mean_severe_ve, multiplier_adj),
              age_group_multiplier_adj %>% dplyr::select(age_group, mean_nonsevere_ve, mean_severe_ve),
              with_protection$severe_ve_pred,
              with_protection$nonsevere_ve_pred,
              group_multiplier_adj %>% dplyr::select(age_group, risk_group, mean_nonsevere_ve, mean_severe_ve)))
}

protection_at_model_init <- protection_at_model_initialization(clean_df[[1]]) 
mult_adj <- protection_at_model_init[[1]]

severe_infection_multipliers <- setDT(merge(merge(average_severe_incidence, average_nonsevere_incidence, by = c("age_group", "risk_group"), all.x = TRUE) %>% 
                                              mutate(multiplier = severe_inc/nonsevere_inc),
                                                  mult_adj, by = c("age_group", "risk_group"), all.x = TRUE) %>%
  mutate(multiplier = if_else(age_group == "0-17 years", 0, multiplier)))

setkeyv(severe_infection_multipliers, c("age_group", "risk_group"))
############################################################################################
#Beta calculation: We have 3 age-specific lambdas (0-17 years, 18-64 years, 65+ years). For the 18-64 year and 65+ year
# lambdas, we need to account for differences in magnitude of the nonsevere incidence and (1-PE) estimate among the further-
# stratified age groups. We are assuming that the 18-64 yr lambda has the 18-29 year age group as baseline, and the 65+ yr
# lambda has the 65-74 year age group as baseline. Beta is fixed. We are using the average nonsevere incidence and average 
# (1-PE) estimates by age group at t = 0 to calculate beta.
beta_calc <- merge(protection_at_model_init[[5]], average_nonsevere_incidence, by = c("age_group", "risk_group"), all.x = TRUE) 

beta_0_17 <- beta_calc %>% filter(age_group == "0-17 years") %>%
  mutate(pe_adj = (1 - mean_nonsevere_ve[1])/(1 - mean_nonsevere_ve),
         inc_adj = nonsevere_inc/nonsevere_inc[1],
         beta = pe_adj * inc_adj) %>%
  dplyr::select(age_group, risk_group, beta, mean_nonsevere_ve)
beta_18_64 <- beta_calc %>% filter(age_group %in% c("18-29 years", "30-49 years", "50-64 years")) %>%
  mutate(pe_adj = (1 - mean_nonsevere_ve[1])/(1 - mean_nonsevere_ve),
         inc_adj = nonsevere_inc/nonsevere_inc[1],
         beta = pe_adj * inc_adj) %>%
  dplyr::select(age_group, risk_group, beta, mean_nonsevere_ve)
beta_65_plus <- beta_calc %>% filter(age_group %in% c("65-74 years", "75+ years")) %>%
  mutate(pe_adj = (1 - mean_nonsevere_ve[1])/(1 - mean_nonsevere_ve),
         inc_adj = nonsevere_inc/nonsevere_inc[1],
         beta = pe_adj * inc_adj) %>%
  dplyr::select(age_group, risk_group, beta, mean_nonsevere_ve)

beta <- setDT(rbind(beta_0_17, beta_18_64, beta_65_plus))
setkeyv(beta, c("age_group", "risk_group"))


calibration <- function() {
  
  total_inf <- sum(inf_by_age$total_inf) 
  
  avg_nonsevere_inc_adj <- average_nonsevere_incidence %>%
    mutate(nonsevere_inc = (nonsevere_inc/100000)/7)
  
  lambda_0_17 <- sum((avg_nonsevere_inc_adj %>% filter(age_group == "0-17 years"))$nonsevere_inc)/sum((beta%>% filter(age_group == "0-17 years"))$beta  * (1-(beta%>% filter(age_group == "0-17 years"))$mean_nonsevere_ve) * (total_inf/50000000))
  lambda_18_64 <- sum((avg_nonsevere_inc_adj %>% filter(age_group %in% c("18-29 years", "30-49 years", "50-64 years")))$nonsevere_inc)/
    sum((beta%>% filter(age_group %in% c("18-29 years", "30-49 years", "50-64 years")))$beta  * (1-(beta%>% filter(age_group %in% c("18-29 years", "30-49 years", "50-64 years")))$mean_nonsevere_ve) * (total_inf/50000000))
  lambda_65 <- sum((avg_nonsevere_inc_adj %>% filter(age_group %in% c("65-74 years", "75+ years")))$nonsevere_inc)/sum((beta%>% filter(age_group %in% c("65-74 years", "75+ years")))$beta  * (1-(beta%>% filter(age_group %in% c("65-74 years", "75+ years")))$mean_nonsevere_ve) * (total_inf/50000000))
  
  return(c(lambda_0_17, lambda_18_64, lambda_65))
}

lambdas <- calibration()
