########################################################################################################################
#Title: Vaccination Model
#Author: Hailey Park
#Date: February 5th, 2024
########################################################################################################################

# Function to calculate infection risk and simulate outcomes for all individuals
# Returns: list of 4 vectors (length = pop_size each)
#   1. severe_outcomes: binary indicator of severe infection
#   2. nonsevere_outcomes: binary indicator of non-severe infection  
#   3. severe_pe: updated protection against severe outcomes (0-1)
#   4. nonsevere_pe: updated protection against non-severe outcomes (0-1)
outcome_occurrence <- function(age, inf, day, risk, immuno, vacc, perfect_immunity_counter, inf_tracker_df, lambda, betas_vec, 
                               severe_mult_vec, prior_protection_severe, prior_protection_nonsevere, prior_time_since, individual) {
  
  # Initialize protection levels with prior values
  severe_pe <- prior_protection_severe
  nonsevere_pe <- prior_protection_nonsevere
  
  # ============================================================================
  # Calculate dynamic contact term for force of infection (i.e., C_ij * I_i/N_i)
  # ============================================================================
  
  # Contact matrix weighted by: (infections in past 5 days / population) * contact rate
  # Calculate the age contact matrix terms for all individuals (dynamic term)
  #     If no circulating infections, return no risk
  # Note: Uses rows 4-8 of inf_tracker_df to account for 3-day latent period
  contact_matrix_term_by_age <- data.table(
    age_group = c("0-17 years","18-29 years", "30-49 years", "50-64 years", "65-74 years", "75+ years"),
    contact_matrix_term = c(
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X0.17.years),
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X18.29.years),
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X30.49.years),
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X50.64.years),
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X65.74.years),
      sum(colSums(inf_tracker_df[4:8, -1])/inf_by_age$total_pop * contact_matrix$X75..years)
    )
  )
  
  # Early exit: if no circulating infections, return no outcomes
  if (sum(contact_matrix_term_by_age$contact_matrix_term) == 0) {
    return(list(rep(0, length(age)), rep(0, length(age)), 
                rep(0, length(age)), rep(0, length(age))))
  }
  
  # Map contact matrix terms to each individual based on their age group
  daily_contact_matrix_term <- data.table(individual = individual, age_group = age)[
    contact_matrix_term_by_age,
    on = "age_group",
    nomatch = NULL
  ] %>%
    arrange(individual) %>%
    pull(contact_matrix_term)
  
  # ============================================================================
  # Update individual's current protection levels against severe & non-severe infection
  #     Stored in vectors severe_pe and non_severe_pe
  # ============================================================================
  
  #Identify individuals eligible for infection and whose waning curves need to be updated to merge with 
  # waning_data_clean to get protection at specific time point
  # NOTE: Eligible individuals are those who do not active perfect immunity, and those whose immunity hasn't been updated in 14 days 
  #         (i.e., who's prior time-since is different from current time-since under the function (ceiling(X/14)))
  index_individuals_eligible <- which(
    perfect_immunity_counter == 0 & 
      (ceiling((prior_time_since)/14) != ceiling(day/14))
    )
  
  # Create lookup table for eligible individuals to merge with waning curves
  df_individuals_eligible <- data.table(
    index_individual = index_individuals_eligible,
    age_group = age[index_individuals_eligible],
    risk_group = risk[index_individuals_eligible],
    prior_inf = if_else(inf[index_individuals_eligible] == "noinf", 0, 1),
    immuno = immuno[index_individuals_eligible],
    prior_vacc = if_else(vacc[index_individuals_eligible] == "unvax", 0, 1),
    days = day[index_individuals_eligible], # time_since_last
    weeks = pmin(ceiling(day[index_individuals_eligible]/14), 52), # Converting day into week
    key = c("immuno", "weeks", "prior_inf", "prior_vacc")
  )[order(index_individual, decreasing=FALSE), ]
  
  # # KMB - I don't think we need these lines because severe_pe is already a copy of prior_protection_severe.. check later
  # # For individuals whose time-since is still in same 2-week bin as prior time-since, 
  # # keep prior protection (i.e., only update every 14 days)
  # index_individuals_no_updated_waning <- which(
  #   perfect_immunity_counter == 0 & 
  #     (ceiling((prior_time_since)/14) == ceiling(day/14))
  # )
  # severe_pe[index_individuals_no_updated_waning] <- prior_protection_severe[index_individuals_no_updated_waning]
  # nonsevere_pe[index_individuals_no_updated_waning] <- prior_protection_nonsevere[index_individuals_no_updated_waning]
  
  # Merge eligible individuals with waning curve data to get updated protection
  df_individuals_eligible <- merge(df_individuals_eligible, waning_data_clean_alt, all.x = TRUE)
  setkeyv(df_individuals_eligible, c("age_group", "risk_group"))
  
  #print(df_individuals_eligible %>% filter(index_individual == 2140979)) #print out individual waning protection throughout simulation as a check
  #print(df_individuals_eligible %>% filter(index_individual == 1325471))
  
  # Immune-naive individuals (no vaccination, no prior infection) have no protection
  immune_naive_index <- which(
    df_individuals_eligible$prior_vacc == 0 & 
      df_individuals_eligible$prior_inf == 0
  )
  df_individuals_eligible[immune_naive_index, c("severe_ve", "nonsevere_ve")] <- 0
  
  # Update protection levels (severe & non-severe) for eligible individuals
  severe_pe[df_individuals_eligible$index_individual] <- df_individuals_eligible$severe_ve
  nonsevere_pe[df_individuals_eligible$index_individual] <- df_individuals_eligible$nonsevere_ve
  
  # Set protection for individuals in perfect immunity window (90 days post-infection)
  severe_pe[which(perfect_immunity_counter > 0)] <- 1
  nonsevere_pe[which(perfect_immunity_counter > 0)] <- 1
  
  # ============================================================================
  # Calculate infection risk and simulate outcomes
  # ============================================================================
  
  # Calculate risk
  nonsevere_risk <- lambda * betas_vec * (1 - nonsevere_pe) * daily_contact_matrix_term
  severe_risk <- lambda * betas_vec * (1 - severe_pe) * daily_contact_matrix_term * severe_mult_vec
  
  # Convert to probability P(infection) = 1 - exp(risk)
  nonsevere_prob <- 1 - exp(-nonsevere_risk)
  severe_prob <- 1 - exp(-severe_risk)
  
  # Simulate outcomes
  severe_outcomes <- rbinom(length(severe_prob), 1, severe_prob)
  nonsevere_outcomes <- rbinom(length(nonsevere_prob), 1, nonsevere_prob)
  
  # Return list of outcome vectors and updated protection levels
  return(list(severe_outcomes, nonsevere_outcomes, severe_pe, nonsevere_pe))
}

########################################################################################################################
#These are the simulations for different vaccine frequencies (annual, semi-annual) over an 18-month period (July 1, 2023 - January 1, 2025)
# NOTE: Even though the pre-assignment of vaccine uptake under different vaccine scenarios already accounts for simulating different 
#       frequencies of vaccination (annual, semi-annual), the reason I wrote these different functions is because we are also estimating annual 
#       risk among those vaccinated, only if they receive the full number of doses they are eligible for. 
#       - For the annual strategies, it's simply counting outcomes among individuals with both `vax_1` and `vax_3`,
#       - For semi-annual strategies 5 + 7, it's simply counting outcomes among individuals with vax_1`, `vax_2`, and `vax_3`
#       - For semi-annual strategies 8-10 + historical, some individuals only eligible for `vax_1` and `vax_3`, others are eligible 
#         for all 3. So must take that into account

simulation_semiannual_strat_8_9_10 <- function(input_df, sim_days, rel_ve_check = FALSE) {
  
  # ============================================================================
  # Initialization
  # ============================================================================
  
  # Initialize data frame to store severe and non-severe outcome counts, stratified by age and risk group 
  #     NOTE: The data frame is wide because it is storing outcome counts for severe and non-severe infections separately across 665 days (18-month simulation)
  grouped_outcome_counts <- input_df  %>% 
    group_by(age_group, risk_group) %>% 
    summarise(total_pop = n())
  grouped_outcome_counts[sprintf("day%s", (1:665))] <- NA
  grouped_outcome_counts[sprintf("nonsevere_day%s", (1:665))] <- NA
  
  # Store protective effectiveness (for debugging)
  grouped_PE_counts <- input_df  %>% 
    group_by(age_group, risk_group) %>% 
    summarise(total_pop = n())
  grouped_PE_counts[sprintf("day%s",(1:665))] <- NA
  grouped_PE_counts[sprintf("nonsevere_day%s",(1:665))] <- NA
  
  #Input data (entire pop)
  input_df <- input_df %>% arrange(individual)
  input_df[sprintf("day%s", (1))] <- NA
  input_df[sprintf("nonsevere_day%s", (1))] <- NA
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each time step
  individual <- input_df$individual
  age <- as.character(input_df$age_group)
  risk <- as.character(input_df$risk_group)
  vacc <- as.character(input_df$prior_vacc)
  inf <- as.character(input_df$prior_inf)
  immuno <- input_df$immuno
  time_since_last <- pmax(input_df$days_since_last_dose_inf, 1)
  time_since_last_dose <- pmax(input_df$days_since_last_dose, 1)
  
  #Initialize counter to track perfect immunity for each individual (vector pop_size x 1)
  #   Individuals infected in 3 months preceding start of sim have perfect immunity at start
  #   Assumed length of perfect immunity is 90 days
  perfect_immunity_counter <- rep(0, nrow(input_df))
  index_recent_infection <- which(inf != "noinf" & 
                                    time_since_last < 90 & 
                                    ((time_since_last < time_since_last_dose) | (is.na(time_since_last_dose)))) 
  perfect_immunity_counter[index_recent_infection] <- 91 - time_since_last[index_recent_infection] 
  
  #Initialize daily age-specific infection trackers using initial conditions stored in inf_by_age
  daily_infection_by_age <- inf_by_age$total_inf
  inf_tracker_df <- infection_tracker_df
  
  #Initialize vectors for vaccine waves (vectors pop_size x 1)
  vaccine_wave1 <- input_df$vaccine_wave1
  vaccine_wave2 <- input_df$vaccine_wave2
  vaccine_wave3 <- input_df$vaccine_wave3
  vaccine_wave4 <- input_df$vaccine_wave4
  
  #Assign individual beta 
  #   These are age-specific adjustments based on differences in magnitude of non-severe incidence and (1 - PE_non-severe)), 
  #   and lambda calculated at t = 0
  betas_vec <- data.table(individual = individual, age_group = age)[
    beta,
    on = "age_group",
    nomatch = NULL
  ] %>%
    arrange(individual) %>%
    pull(betas)
  
  #Assign age- and risk-specific severe multipliers
  #   This multiplier is the inverse of non-severe to severe incidence
  severe_multiplier_vec <- data.table(
    individual = individual, 
    age_group = age, 
    risk_group = risk
  )[
    severe_infection_multipliers,
    on = c("age_group", "risk_group"),
    nomatch = NULL
  ] %>%
    arrange(individual) %>%
    pull(multiplier)
  
  #Initialize vector to store prior time-since and protection, and set it to values at model initialization
  #   NOTE: We update individual's immunity every 14 days instead of every day to improve run time
  prior_time_since <- time_since_last - 1
  prior_time_since[prior_time_since >= 730] <- 730     #Assuming that >24 month waning is same as 24 month waning PE
  prior_time_since[prior_time_since < 0] <- 0
  prior_protection_severe <- protection_at_model_init[[2]] 
  prior_protection_nonsevere <- protection_at_model_init[[3]] 
  
  #Initialize vectors to store vax outcomes by age and risk group
  vax_outcomes_severe <- rep(0, 18)
  vax_outcomes_nonsevere <- rep(0, 18)
  
  # For debugging - KMB
  mean_protection <- data.frame(matrix(ncol = 2, nrow = 92)) 
  
  # Define month boundaries (day numbers from July 1, 2023)
  #   Day 1 = July 1, 2023
  
  #   Days 1-31: July 2023 - lambda1
  #   Days 32-62: August 2023 - lambda2
  #   Days 63-92: September 2023 - lambda3
  #   Days 93-123: October 2023 - lambda4
  #   Days 124-153: November 2023 - lambda5
  #   Days 154-184: December 2023 - lambda6
  #   Days 185-215: January 2024 - lambda7
  #   Days 216-243: February 2024 - lambda8
  #   Days 244-274: March 2024 - lambda9
  #   Days 275-304: April 2024 - lambda10
  #   Days 305-335: May 2024 - lambda11
  #   Days 336-365: June 2024 - lambda12
  #   Days 366-396: July 2024 - lambda13
  #   Days 397-427: August 2024 - lambda14
  #   Days 428-457: September 2024 - lambda15
  #   Days 458-488: October 2024 - lambda16
  #   Days 489-518: November 2024 - lambda17
  #   Days 519-549: December 2024 - lambda18
  #   Days 550-580: January 2025 - lambda19
  #   Days 581-608: February 2025 - lambda20
  #   Days 609-639: March 2025 - lambda21
  #   Days 640-665: April 2025 - lambda22
  
  month_boundaries <- c(1, 32, 63, 93, 124, 154, 185, 216, 244, 275, 305, 336, 
                        366, 397, 428, 458, 489, 519, 550, 581, 609, 640, 666)
  
  lambda_values <- list(lambda_1, lambda_2, lambda_3, lambda_4, lambda_5, lambda_6,
                        lambda_7, lambda_8, lambda_9, lambda_10, lambda_11, lambda_12,
                        lambda_13, lambda_14, lambda_15, lambda_16, lambda_17, lambda_18,
                        lambda_19, lambda_20, lambda_21, lambda_22)
  
  # FIXME: update for vax_4
  # KMB - double check this section AND I think rough estimate of rel. VE calculation will come into play here
  # Store vaccination status (fully vaccinated = 1, otherwise = 0)
  #     For ages 65-74 and 75+: requires vax_1, vax_2, and vax_3
  #     For younger ages + healthy/higher risk: requires vax_1 and vax_3 only
  input_df$vax_status <- as.integer(as.logical(input_df$vax_1 == 1 & 
                                                 input_df$vax_2 == 1 & 
                                                 input_df$vax_3 == 1))
  input_df$vax_status[which(age %in% c("0-17 years", "18-29 years", "30-49 years", "50-64 years") & 
                              risk %in% c("healthy", "higher risk") & 
                              input_df$vax_1 == 1 & input_df$vax_3 == 1)] <- 1
  
  if (rel_ve_check == TRUE){
    # Store number of severe infections in those with vax 1 and those without vax 1 from March 1 2024-August 31 2024
    # Days 244 - 427
    
    #vax1_status <- as.integer(as.logical(vacc == "boosted_2nd" & input_df$vax_1 == 1))
    #novax1_status <- as.integer(as.logical(vacc == "boosted_2nd" & input_df$vax_1 == 0))
    
    rel_ve_df <- data.frame(severe_vax1 = 0,
                            severe_no_vax1 = 0,
                            #total_vax1 = sum(vax1_status),
                            #total_no_vax1 = sum(novax1_status),
                            total_vax1 = 0,
                            total_no_vax1 = 0)
  }
  
  #sum(vax1_status)
  #sum(boosted_novax1_status)
  # if there's enough people who received primary series but not this most recent booster, compare attack rate in that grp instead
  
  
  # ============================================================================
  # Iterate through each time step
  # ============================================================================
  for (i in (1:sim_days)) {
    
    print(paste0("Day: ", i)) 
    
    # ============================================================================
    # Update time_since_last for those who receive a vaccine on day i
    #     Vaccine timing is stored in vaccine_wave vectors
    #     New protection term will be computed later in outcome_occurrence function
    # ============================================================================
    
    # Vax_1 staggered over first 365 days
    if (i %in% c(1:365)){
      vaccine_wave1_index <- which(vaccine_wave1 == i)
      time_since_last[vaccine_wave1_index] <- 1
    }
    
    # Vax_2 staggered over the last 126 days of year 1
    if (i %in% c(240:365)){
      vaccine_wave2_index <- which(vaccine_wave2 == i)
      time_since_last[vaccine_wave2_index] <- 1
    }
    
    # Vax_3 staggered over 182 days in fall 2024
    if (i %in% c((366:547))){
      vaccine_wave3_index <- which(vaccine_wave3 == i)
      time_since_last[vaccine_wave3_index] <- 1
    }
    
    # Vax_4 staggered spring 2025
    if (i %in% c((609:665))){
      vaccine_wave4_index <- which(vaccine_wave4 == i)
      time_since_last[vaccine_wave4_index] <- 1
    }
    
    #Assuming that >24 month waning is same as 24 month waning PE
    time_since_last[time_since_last >= 730] <- 730     
    time_since_last[time_since_last <= 0] <- 1
    
    # ============================================================================
    # Pull month-specific lambda (depending on what month i falls into)
    # ============================================================================
    
    month_index <- findInterval(i, month_boundaries)
    lambdas <- lambda_values[[month_index]]
    
    # ============================================================================
    # Simulate infections on day i
    # ============================================================================
    # Generate outcomes for all individuals
    outcomes <- outcome_occurrence(age, inf, time_since_last, risk, immuno, vacc, perfect_immunity_counter, 
                                   inf_tracker_df, lambdas, betas_vec, severe_multiplier_vec, 
                                   prior_protection_severe, prior_protection_nonsevere, prior_time_since, individual)
    
    # Extract outcome vectors (each length = pop_size)
    severe_outcomes             <- outcomes[[1]] # Binary: severe infection indicator
    nonsevere_outcomes          <- outcomes[[2]] # Binary: non-severe infection indicator
    prior_protection_severe     <- outcomes[[3]] # Numeric [0-1]: updated protection against severe
    prior_protection_nonsevere  <- outcomes[[4]] # Numeric [0-1]: updated protection against non-severe
    
    print(paste0("Total daily severe infections: ", sum(severe_outcomes)))
    print(paste0("Total daily nonsevere infections: ", sum(nonsevere_outcomes)))
    
    # ============================================================================
    # Update infection status and time counters
    # ============================================================================
    
    # Store previous time-since values before updating (vector length = pop_size)
    prior_time_since <- time_since_last
    
    # For individuals with no infection today, increment days since last infection
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1
    
    # Decrement perfect immunity counter by 1 day (for those still in 90-day window)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1
    
    # For individuals with new infections today:
    # - Mark as previously infected (used in look up table in outcome occurrence when assigning corresponding waning curve)
    # - Reset time since last to 1 day
    # - Reset perfect immunity to 90 days
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1 
    time_since_last[index_outcome] <- 1
    perfect_immunity_counter[index_outcome] <- 90
    
    # Handle edge case: if individual has both severe and non-severe infection, keep only severe
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0
    
    # ============================================================================
    # Update infection tracking by age group
    # ============================================================================
    # Count total infections (severe + non-severe) by age group for today
    daily_infection_by_age[1] <- sum(severe_outcomes[which(age == "0-17 years")]) + sum(nonsevere_outcomes[which(age == "0-17 years")])
    daily_infection_by_age[2] <- sum(severe_outcomes[which(age == "18-29 years")]) + sum(nonsevere_outcomes[which(age == "18-29 years")])
    daily_infection_by_age[3] <- sum(severe_outcomes[which(age == "30-49 years")]) + sum(nonsevere_outcomes[which(age == "30-49 years")])
    daily_infection_by_age[4] <- sum(severe_outcomes[which(age == "50-64 years")]) + sum(nonsevere_outcomes[which(age == "50-64 years")])
    daily_infection_by_age[5] <- sum(severe_outcomes[which(age == "65-74 years")]) + sum(nonsevere_outcomes[which(age == "65-74 years")])
    daily_infection_by_age[6] <- sum(severe_outcomes[which(age == "75+ years")]) + sum(nonsevere_outcomes[which(age == "75+ years")])
    
    # Update infection tracker: rolling 8-day window
    # inf_tracker_df: 8 rows x (1 + 6 age groups), tracks infections in past 8 days
    inf_tracker_df <- inf_tracker_df[1:7, ]
    inf_tracker_df <- rbind(c(1, daily_infection_by_age), inf_tracker_df)
    inf_tracker_df$days_since <- c(1:8)
    
    # ============================================================================
    # Record outcomes in output data frame
    # ============================================================================
    # Add today's individual-level outcomes to data frame
    input_df$thisday <- severe_outcomes
    input_df$nonsevere_thisday <- nonsevere_outcomes
    
    # Count infections among vaccinated individuals
    input_df$vax_thisday <- as.integer(as.logical(severe_outcomes == 1 & input_df$vax_status == 1))
    input_df$vax_nonsevere_thisday <- as.integer(as.logical(nonsevere_outcomes == 1 & input_df$vax_status == 1))
    
    # ============================================================================
    # Aggregate outcomes by age and risk group
    # ============================================================================
    
    # Summarize outcomes by age_group x risk_group strata
    grouped_outcomes_today <- input_df %>% 
      group_by(age_group, risk_group) %>% 
      summarise(
        total_severe = sum(thisday),
        total_nonsevere = sum(nonsevere_thisday),
        total_severe_vax = sum(vax_thisday),
        total_nonsevere_vax = sum(vax_nonsevere_thisday)
      )
    
    # Store daily counts in output matrix
    # Columns 4 to 668: severe outcomes by day (i=1 to i=665)
    # Columns 669 to FIXME: non-severe outcomes by day
    grouped_outcome_counts[, i + 3] <- grouped_outcomes_today$total_severe
    grouped_outcome_counts[, i + (665 + 3)] <- grouped_outcomes_today$total_nonsevere
    
    # Accumulate vaccine breakthrough infections across all days
    vax_outcomes_severe <- vax_outcomes_severe + grouped_outcomes_today$total_severe_vax
    vax_outcomes_nonsevere <- vax_outcomes_nonsevere + grouped_outcomes_today$total_nonsevere_vax
    
    # If checking relative VE
    if (rel_ve_check == TRUE){
      # Days 244 - 427
      # i == 365 means ending in June 2024
      # i == 427 means ending in August 2024
      if (i == 244){
        vax1_status <- as.integer(vacc %in% c("boosted_2nd") & #"fullvax", "boosted_1st", "boosted_2nd") &
                                    risk == "healthy" &
                                    #inf == "no_inf" & 
                                    input_df$vax_1 == 1)
        
        novax1_status <- as.integer(vacc %in% c("boosted_2nd") & #"fullvax", "boosted_1st", "boosted_2nd") &
                                      risk == "healthy" &
                                      #inf == "no_inf" & 
                                      input_df$vax_1 == 0)
      }
      
      if (i > 243) {
        rel_ve_df$severe_vax1 <- rel_ve_df$severe_vax1 + sum(as.integer(as.logical(severe_outcomes == 1 & vax1_status == 1)))
        rel_ve_df$severe_no_vax1 <- rel_ve_df$severe_no_vax1 + sum(as.integer(as.logical(severe_outcomes == 1 & novax1_status == 1)))
      }
      
      if (i == 427){
        rel_ve_df$total_vax1 = sum(vax1_status)
        rel_ve_df$total_no_vax1 = sum(novax1_status)
        
        return(rel_ve_df)
      }
    }
    
    # On final day (day 665), record cumulative vaccination statistics
    if (i == 665) {
      grouped_outcomes_today <- input_df %>% 
        group_by(age_group, risk_group) %>% 
        summarise(
          total_vaccines = sum(vax_1) + sum(vax_2) + sum(vax_3),
          total_vaccinated = sum(vax_status)
        )
      
      grouped_outcome_counts$total_vaccines <- grouped_outcomes_today$total_vaccines
      grouped_outcome_counts$total_vaccinated <- grouped_outcomes_today$total_vaccinated
      grouped_outcome_counts$total_severe_vax <- vax_outcomes_severe
      grouped_outcome_counts$total_nonsevere_vax <- vax_outcomes_nonsevere
    }
    
    # For debugging -- KMB
    input_df$PE_severe <- prior_protection_severe
    input_df$PE_nonsevere <- prior_protection_nonsevere
    
    grouped_PE <- input_df %>%
      group_by(age_group, risk_group) %>%
      summarise(
        severePE = mean(PE_severe),
        nonseverePE = mean(PE_nonsevere)
      )
    
    grouped_PE_counts[, i + 3] <- grouped_PE$severePE
    grouped_PE_counts[, i + (665 + 3)] <- grouped_PE$nonseverePE
    
    mean_protection[i,] <- c(mean(prior_protection_severe), mean(prior_protection_nonsevere)) # KMB
  }
  
  return(grouped_outcome_counts)
  #return(grouped_PE_counts)
}
