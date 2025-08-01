########################################################################################################################
#Title: Vaccination Interventions (Dynamic with Age Mixing)
#Author: Hailey Park
#Date: February 5th, 2024
########################################################################################################################

#Function for outcome occurrence based on individual-level risk (Risk = Lambda * (1 - PE))
outcome_occurrence <- function(age, inf, day, risk, immuno, vacc, perfect_immunity_counter, inf_tracker_df, contact_matrix_adj_factors, lambda, betas, severe_mult, prior_protection_severe, prior_protection_nonsevere, prior_time_since, individual) {
  
  severe_pe <- prior_protection_severe
  nonsevere_pe <- prior_protection_nonsevere
  
  #Calculate the age contact matrix terms for all individuals (dynamic term) -- if no circulating infections, return no risk (optimizing function)
  contact_matrix_term_by_age <- data.table(age_group = c("0-17 years","18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),
                                           contact_matrix_term = c(sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X0.17.years),
                                                                   sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X18.29.years),
                                                                   sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X30.49.years),
                                                                   sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X50.64.years),
                                                                   sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X65.74.years),
                                                                   sum(colSums(inf_tracker_df[4:8,-1])/inf_by_age$total_pop * contact_matrix$X75..years)))

    if (sum(contact_matrix_term_by_age$contact_matrix_term) == 0){
    return(list(rep(0, length(age)), rep(0, length(age)), rep(0, length(age)), rep(0, length(age))))
  }
  
  #Assign contact matrix term to each individual
  daily_contact_matrix_term <- ((data.table(individual = individual, age_group = age)[contact_matrix_term_by_age, 
                                                            on=c("age_group"), 
                                                            nomatch = NULL]) %>% arrange(individual))$contact_matrix_term
  
  #Creating a dataframe of individuals eligible for infection to merge with waning_data_clean to get protection at specific time point
  # NOTE: Eligible individuals are those who do not active perfect immunity, and who's prior time-since is different from current time-since
  #       under the function (ceiling(X/14))
  index_individuals_eligible <- which(perfect_immunity_counter == 0 & (ceiling((prior_time_since)/14) != ceiling(day/14)))
  df_individuals_eligible <- data.table(index_individual = index_individuals_eligible,
                                        age_group = age[index_individuals_eligible],
                                        risk_group = risk[index_individuals_eligible],
                                        prior_inf = if_else(inf[index_individuals_eligible] == "noinf", 0, 1),
                                        immuno = immuno[index_individuals_eligible],
                                        prior_vacc = if_else(vacc[index_individuals_eligible] == "unvax", 0, 1),
                                        days = day[index_individuals_eligible],
                                        weeks = pmin(ceiling(day[index_individuals_eligible]/14), 52), #converting day into week
                                        key = c("immuno", "weeks", "prior_inf", "prior_vacc"))[order(index_individual,decreasing=FALSE),]
  
  #For individuals whose time-since are within same bi-week as prior time-since, no need to update waning estimate
  index_individuals_no_updated_waning <- which(perfect_immunity_counter == 0 & (ceiling((prior_time_since)/14) == ceiling(day/14)))
  severe_pe[index_individuals_no_updated_waning] <- prior_protection_severe[index_individuals_no_updated_waning]
  nonsevere_pe[index_individuals_no_updated_waning] <- prior_protection_nonsevere[index_individuals_no_updated_waning]

  #For individuals who time-since are not within the same bi-week as prior time-since, update waning estimate
  df_individuals_eligible <- merge(df_individuals_eligible, waning_data_clean_alt, all.x = TRUE)
  setkeyv(df_individuals_eligible, c("age_group", "risk_group"))

  # print(df_individuals_eligible %>% filter(index_individual == 2140979))
  #print(df_individuals_eligible %>% filter(index_individual == 1325471))

  #Individuals who are unvaccinated and no prior infection history has no protection
  immune_naive_index <- which(df_individuals_eligible$prior_vacc == 0 & df_individuals_eligible$prior_inf == 0)
  df_individuals_eligible[immune_naive_index, c("severe_ve", "nonsevere_ve")] <- 0

  #Calculate protection (severe + nonsevere) for individuals eligible for infection and updating waning
  severe_pe[df_individuals_eligible$index_individual] <- df_individuals_eligible$severe_ve
  nonsevere_pe[df_individuals_eligible$index_individual] <- df_individuals_eligible$nonsevere_ve

  #Individuals with perfect immunity, set immunity to 100%
  severe_pe[which(perfect_immunity_counter > 0)] <- 1
  nonsevere_pe[which(perfect_immunity_counter > 0)] <- 1
  
  #Calculate risk
  nonsevere_risk <- lambda * betas * (1 - nonsevere_pe) * contact_matrix_adj_factors * (daily_contact_matrix_term)
  severe_risk <- lambda * betas * (1 - severe_pe)  *  contact_matrix_adj_factors * (daily_contact_matrix_term) * severe_mult
  
  #if nonsevere or severe risk > 1, set to 1
  nonsevere_risk[nonsevere_risk > 1] <- 1
  severe_risk[severe_risk > 1] <- 1
  
  #Simulate outcomes
  severe_outcomes <- rbinom(length(severe_risk), 1, severe_risk)
  nonsevere_outcomes <- rbinom(length(nonsevere_risk), 1, nonsevere_risk)
  
  return(list(severe_outcomes, nonsevere_outcomes, severe_pe, nonsevere_pe))
}
########################################################################################################################
#This is the one year booster simulation, July 2023-July 2024

historicalVaccinationSimulation <- function(df, params){
  
  #Set lambda values from `params` vector
  lambda_1 <- params[1]
  lambda_2 <- params[2]
  lambda_3 <- params[3]
  lambda_4 <- params[4]
  lambda_5 <- params[5]
  lambda_6 <- params[6]
  lambda_7 <- params[7]
  lambda_8 <- params[8]
  lambda_9 <- params[9]
  lambda_10 <- params[10]
  lambda_11 <- params[11]
  lambda_12 <- params[12]
  lambda_13 <- params[13]
  lambda_14 <- params[14]
  lambda_15 <- params[15]
  lambda_16 <- params[16]
  lambda_17 <- params[17]
  lambda_18 <- params[18]
  
  
  lambda_0_17 <- params[21]
  lambda_18_64 <- params[22]
  lambda_65_plus <- params[23]
  
  #Store severe and nonsevere outcome counts in grouped dataframe, stratified by age and risk group 
  # NOTE: The dataframe is wide because it is storing outcome counts for severe and non-severe infections separately across 547 days (18-month simulation)
  grouped_outcome_counts <- df  %>% group_by(age_group, risk_group) %>% summarise(total_pop = n())
  grouped_outcome_counts[sprintf("day%s",(1:547))] <- NA
  grouped_outcome_counts[sprintf("nonsevere_day%s",(1:547))] <- NA
  
  #Input data (entire pop)
  input <- df %>% arrange(individual)
  
  #Set empty columns to populate for severe/nonsevere outcomes
  input$day1 <- NA
  input$nonsevere_day1 <- NA
  
  #Population's info (age_group, num_doses, prior_inf, etc.)
  individual <- input$individual
  age <- as.character(input$age_group)
  risk <- as.character(input$risk_group)
  vacc <- as.character(input$prior_vacc)
  inf <- as.character(input$prior_inf)
  immuno <- input$immuno
  time_since_last <- pmax(input$days_since_last_dose_inf, 1)
  time_since_last_dose <- pmax(input$days_since_last_dose, 1)
  
  #If infection occurs, counting down perfect immunity (90 days)
  perfect_immunity_counter <- rep(0,nrow(input))
  
  #Individuals infected in 3 months preceding start of sim have perfect immunity at start
  index_recent_infection <- which(inf != "noinf" & time_since_last < 90 & ((time_since_last < time_since_last_dose) | (is.na(time_since_last_dose)))) 
  perfect_immunity_counter[index_recent_infection] <- 91 - time_since_last[index_recent_infection] 
  
  #Set daily age-specific infection trackers to values from model initialization
  daily_infection_by_age <- inf_by_age$total_inf
  inf_tracker_df <- infection_tracker_df
  
  #Set vectors for vaccine waves
  vaccine_wave <- input$vaccine_wave
  second_vaccine_wave <- input$second_vaccine_wave
  third_vaccine_wave <- input$third_vaccine_wave
  
    
  #Assign age-specific lambda to each individual
  age_specific_lambdas <- (data.table(individual = individual, age_group = age) %>%
                             mutate(lambda = case_when(age == "0-17 years" ~ lambda_0_17,
                                    age %in% c("18-29 years", "30-49 years", "50-64 years") ~ lambda_18_64,
                                    age %in% c("65-74 years", "75+ years") ~ lambda_65_plus,
                                    TRUE ~ NA)) %>% arrange(individual))$lambda
  
  #Assign beta (These are age-specific lambda adjustments to the age-specific lambdas based on differences in magnitude of nonsevere incidence and (1-PE))
  betas <- ((data.table(individual = individual, age_group = age, risk_group = risk)[beta,
                                        on=c("age_group", "risk_group"),
                                        nomatch = NULL]) %>% arrange(individual))$beta
  

  #Assign severe multipliers (This multiplier is the inverse of nonsevere to severe incidence, and also accounts for the difference in magnitude between severe and nonsevere protection)
  severe_multiplier_with_adj <- ((data.table(individual = individual, age_group = age, risk_group = risk)[severe_infection_multipliers, 
                                                                                 on=c("age_group", "risk_group"), 
                                                                                 nomatch = NULL]) %>% 
                                   mutate(severe_mult = multiplier/multiplier_adj) %>% arrange(individual))$severe_mult
  
  
  #Contact matrix adjustment @ t = 0 (These are age-specific adjustments to account for differences in magnitude from including social mixing)
  contact_matrix_adj_factors <-((data.table(individual = individual, age_group = age)[contact_matrix_adj,
                                                            on=c("age_group"), 
                                                            nomatch = NULL]) %>% arrange(individual))$contact_matrix_adj
  
  #Creating vector to storing prior time-since and protection, and setting it to values at model initialization
  # NOTE: We do this to optimize runtime by updating individual's immunity every 14 days instead of every day.
  prior_time_since <- time_since_last - 1
  prior_time_since[prior_time_since >= 730] <- 730     #Assuming that >24 month waning is same as 24 month waning pe
  prior_time_since[prior_time_since < 0] <- 0
  prior_protection_severe <- protection_at_model_init[[2]] 
  prior_protection_nonsevere <- protection_at_model_init[[3]] 
  
    #Iterate through each time step
  for (i in (1:547)) {
    
    print(paste0("Day: ", i)) 
    
    #Staggering updated booster vaccination over 365 days
    if(i %in% c(1:365)){
      vaccine_wave_index <- which(vaccine_wave == i)
      time_since_last[vaccine_wave_index] <- 1
    }

    #Staggering updated booster (2nd dose) vaccination between the last 126 days
    if(i %in% c(240:365)){
      vaccine_wave_2_index <- which(second_vaccine_wave == i)
      time_since_last[vaccine_wave_2_index] <- 1
    }
    
    #Staggering updated booster vaccination during model validation period over last 182 days
    if(i %in% c((1 + 365):(365 + 182))){
      vaccine_wave_3_index <- which(third_vaccine_wave == i)
      time_since_last[vaccine_wave_3_index] <- 1
    }
    
    #Modify lambdas based on month
    #Lambda Multiplier 1: (July 1, 2023 - July 31, 2023)
    if(i %in% c(1:31)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_1
    }
    
    #Lambda Multiplier 2: (August 1, 2023 - August 31, 2023)
    if(i %in% c(32:62)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_2
    }
    
    #Lambda Multiplier 3: (September 1, 2023 - September 30, 2023)
    if(i %in% c(63:92)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_3
    }
    
    #Lambda Multiplier 4: (October 1, 2023 - October 31, 2023)
    if(i %in% c(93:123)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_4
    }
    
    #Lambda Multiplier 5: (November 1, 2023 - November 30, 2023)
    if(i %in% c(124:153)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_5
    }
    
    #Lambda Multiplier 6: (December 1, 2023 - December 31, 2023)
    if(i %in% c(154:184)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_6
    }
    
    #Lambda Multiplier 7: (January 1, 2024 - January 31, 2024)
    if(i %in% c(185:215)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_7
    }
    
    #Lambda Multiplier 8: (February 1, 2024 - February 28, 2024)
    if(i %in% c(216:243)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_8
    }
    
    #Lambda Multiplier 9: (March 1, 2024 - March 31, 2024)
    if(i %in% c(244:274)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_9
    }
    
    #Lambda Multiplier 10: (April 1, 2024 - April 30, 2024)
    if(i %in% c(275:304)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_10
    }
    
    #Lambda Multiplier 11: (May 1, 2024 - May 31, 2024)
    if(i %in% c(305:335)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_11
    }
    
    #Lambda Multiplier 12: (June 1, 2024 - June 30, 2024)
    if(i %in% c(336:365)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_12
    }
    
    
    #Lambda Multiplier 13: (July 1, 2024 - July 31, 2024)
    if(i %in% c(366:396)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_13
    }
    
    #Lambda Multiplier 14: (August 1, 2024 - August 31, 2024)
    if(i %in% c(397:427)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_14
    }
    
    #Lambda Multiplier 15: (September 1, 2024 - September 30, 2024)
    if(i %in% c(428:457)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_15
    }
    
    #Lambda Multiplier 16: (October 1, 2024 - October 31, 2024)
    if(i %in% c(458:488)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_16
    }
    
    #Lambda Multiplier 17: (November 1, 2024 - November 30, 2024)
    if(i %in% c(489:518)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_17
    }
    
    #Lambda Multiplier 18: (December 1, 2024 - December 31, 2024)
    if(i %in% c(519:547)) {
      age_specific_lambdas <- age_specific_lambdas * lambda_18
    }
    
    time_since_last[time_since_last >= 730] <- 730     #Assuming that >24 month waning is same as 24 month waning pe
    time_since_last[time_since_last <= 0] <- 1

    #Do outcomes occur?
    outcomes <- outcome_occurrence(age, inf, time_since_last, risk, immuno, vacc, perfect_immunity_counter, inf_tracker_df, contact_matrix_adj_factors, age_specific_lambdas, betas, severe_multiplier_with_adj, prior_protection_severe, prior_protection_nonsevere, prior_time_since, individual)
    severe_outcomes <- outcomes[[1]]
    nonsevere_outcomes <- outcomes[[2]]
    prior_protection_severe <- outcomes[[3]]
    prior_protection_nonsevere <- outcomes[[4]]

    print(paste0("Total daily severe infections: ", sum(severe_outcomes)))
    print(paste0("Total daily nonsevere infections: ", sum(nonsevere_outcomes)))
    
    #set prior_time_since to updated time-since
    prior_time_since <- time_since_last

    #If no outcome occurs, increase time since last
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1

    #Decrease 1 from perfect immunity counter (if applicable)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1

    #If outcome occurs,
    #change their prior infection status to 1, time since last to 1, perfect immunity counter to 90 days
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1
    time_since_last[index_outcome] <- 1
    perfect_immunity_counter[index_outcome] <- 90

    # #Then check if severe outcome is hosp vs. death
    # index_severe_outcome <- which(severe_outcomes == 1)
      
    #If both severe outcome and nonsevere outcome occur in same individual, remove nonsevere outcome
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0

    #Re-update daily_infection_by_age counter with new infection counts
    daily_infection_by_age[1] <- sum(severe_outcomes[which(age == "0-17 years")]) + sum(nonsevere_outcomes[which(age == "0-17 years")])
    daily_infection_by_age[2] <- sum(severe_outcomes[which(age == "18-29 years")]) + sum(nonsevere_outcomes[which(age == "18-29 years")])
    daily_infection_by_age[3] <- sum(severe_outcomes[which(age == "30-49 years")]) + sum(nonsevere_outcomes[which(age == "30-49 years")])
    daily_infection_by_age[4] <- sum(severe_outcomes[which(age == "50-64 years")]) + sum(nonsevere_outcomes[which(age == "50-64 years")])
    daily_infection_by_age[5] <- sum(severe_outcomes[which(age == "65-74 years")]) + sum(nonsevere_outcomes[which(age == "65-74 years")])
    daily_infection_by_age[6] <- sum(severe_outcomes[which(age == "75+ years")]) + sum(nonsevere_outcomes[which(age == "75+ years")])

    #Update inf_tracker_df with new infection counts
    inf_tracker_df <- inf_tracker_df[1:7, ]
    inf_tracker_df <- rbind(c(1, daily_infection_by_age), inf_tracker_df)
    inf_tracker_df$days_since <- c(1:8)
    
    #Add individual-level outcome data to dataframe
    input$day1 <- severe_outcomes
    input$nonsevere_day1 <- nonsevere_outcomes
  
    #Group infection outcomes from the day and store in grouped outcome table
    grouped_outcomes <- input %>% group_by(age_group, risk_group) %>% summarise(total_severe = sum(day1),
                                                                                       total_nonsevere = sum(nonsevere_day1))
    
    grouped_outcome_counts[, i + 3] <- grouped_outcomes$total_severe
    grouped_outcome_counts[, i + (3+547)] <- grouped_outcomes$total_nonsevere
    
  }
  
  return(grouped_outcome_counts)
}
########################################################################################################################
