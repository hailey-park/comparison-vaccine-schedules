########################################################################################################################
# Title: Vaccination assignment under historical coverage
# Author: Hailey Park
# Date: January 16th, 2025
# Last updated: January 28, 2026
# 
# Overview: This function assigns
#       - vaccine doses (vax_1, vax_2, vax_3) 
#       - timing (vaccine_wave1, vaccine_wave2, vaccine_wave3) 
# based on historical age/risk-specific vaccine administration data.
# Semi-annual doses (vax_2) are only assigned to 65+ years and immunocompromised groups.
#
# Inputs:
#   df: data frame number of people x 11 columns that includes prior immune status, age, and risk group
#   first, second, third, and fourth, dose data: data frames with corresponding vaccine coverage by age and risk over time
#
# Output:
#   vax_assignment: same as input df with additional columns for vax_1, vax_2, vax_3 and timing
#########################################################################################################################

historical_vax_assignment <- function(df, dose1_data, dose2_data, dose3_data, dose4_data){
  
  # Initialize vaccination assignment columns
  vax_assignment <- df %>% 
    mutate(vax_1 = 0, vax_2 = 0, vax_3 = 0, vax_4 = 0, 
           vaccine_wave1 = 0, vaccine_wave2 = 0, vaccine_wave3 = 0, vaccine_wave4 = 0)
  
  # ============================================================================
  # Define population indices
  # ============================================================================
  
  # Age group indices
  # Note: For 0-17 years, assume 96.7% are older than 6 months
  age_0_17_index <- which(vax_assignment$age_group == "0-17 years")
  age_0_17_index <- sample(age_0_17_index, length(age_0_17_index) * 0.967)
  
  age_18_29_index <- which(vax_assignment$age_group == "18-29 years")
  age_30_49_index <- which(vax_assignment$age_group == "30-49 years")
  age_50_64_index <- which(vax_assignment$age_group == "50-64 years")
  age_65_74_index <- which(vax_assignment$age_group == "65-74 years")
  age_75_plus_index <- which(vax_assignment$age_group == "75+ years")
  
  # Risk group indices
  risk_healthy_index <- which(vax_assignment$risk_group == "healthy")
  risk_highrisk_index <- which(vax_assignment$risk_group == "higher risk")
  risk_immuno_index <- which(vax_assignment$risk_group == "immunocompromised")
  
  # Previously vaccinated individuals
  prev_vaccinated_index <- which(vax_assignment$prior_vacc != "unvax")
  
  # ============================================================================
  # Assign first vaccine dose (vax_1)
  # ============================================================================
  # Uses age-specific and risk-specific coverage rates
  # Only previously vaccinated individuals are eligible
  
  # Function to assign vaccination (modifies in place, no return for efficiency (because vax_assignment df is large))
  assign_vax_1 <- function(this_age_index, this_risk_index, coverage) {
    eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, prev_vaccinated_index))
    
    if (length(eligible_index) == 0) return(invisible(NULL))
    
    n_target <- length(Reduce(intersect, list(this_age_index, this_risk_index)))
    
    # Assign vaccination probabilistically
    vax_assignment[["vax_1"]][eligible_index] <<- rbinom(
      n = length(eligible_index), 
      size = 1, 
      prob = coverage * n_target / length(eligible_index)
    )
    
    invisible(NULL)
  }
  
  # Age 0-17 years
  assign_vax_1(age_0_17_index, risk_healthy_index, 0.13)
  assign_vax_1(age_0_17_index, risk_highrisk_index, 0.235)
  assign_vax_1(age_0_17_index, risk_immuno_index, 0.235)
  
  # Age 18-29 years
  assign_vax_1(age_18_29_index, risk_healthy_index, 0.0926)
  assign_vax_1(age_18_29_index, risk_highrisk_index, 0.1668)
  assign_vax_1(age_18_29_index, risk_immuno_index, 0.1668)
  
  # Age 30-49 years
  assign_vax_1(age_30_49_index, risk_healthy_index, 0.1308)
  assign_vax_1(age_30_49_index, risk_highrisk_index, 0.2355)
  assign_vax_1(age_30_49_index, risk_immuno_index, 0.2355)

  # Age 50-64 years
  assign_vax_1(age_50_64_index, risk_healthy_index, 0.1439)
  assign_vax_1(age_50_64_index, risk_highrisk_index, 0.2591)
  assign_vax_1(age_50_64_index, risk_immuno_index, 0.2591)
  
  # Age 65-74 years
  assign_vax_1(age_65_74_index, risk_healthy_index, 0.2228)
  assign_vax_1(age_65_74_index, risk_highrisk_index, 0.4010)
  assign_vax_1(age_65_74_index, risk_immuno_index, 0.4010)

  # Age 75+ years
  assign_vax_1(age_75_plus_index, risk_healthy_index, 0.221)
  assign_vax_1(age_75_plus_index, risk_highrisk_index, 0.3979)
  assign_vax_1(age_75_plus_index, risk_immuno_index, 0.3979)

  # ============================================================================
  # Print statements to check assignments of vax 1
  # ============================================================================
  age_groups <- c("0-17 years", "18-29 years", "30-49 years", 
                  "50-64 years", "65-74 years", "75+ years")
  
  targets <- c("0-17 years" = 0.143, 
               "18-29 years" = 0.113, 
               "30-49 years" = 0.159,
               "50-64 years" = 0.217, 
               "65-74 years" = 0.379, 
               "75+ years" = 0.376)
  
  # Checks for overall age group
  cat("\n--- Vaccination rates by age (vax 1) ---\n")
  
  for (age in age_groups) {
    index <- vax_assignment$age_group == age
    vax_rate <- sum(vax_assignment$vax_1[index] == 1) / sum(index)
    
    # Apply adjustment for 0-17 age group
    if (age == "0-17 years") vax_rate <- vax_rate / 0.967
    
    cat(sprintf("Percentage of %s receiving vaccinated: %.3f (target: %.3f)\n", 
                age, vax_rate, targets[age]))
  }
  
  # Checks for age and risk group
  risk_groups <- c("healthy", "higher risk", "immunocompromised")
  
  cat("\n--- Vaccination rates by age and risk group (vax 1) ---\n")
  for (age in age_groups) {
    for (risk in risk_groups) {
      index <- vax_assignment$age_group == age & vax_assignment$risk_group == risk
      vax_rate <- sum(vax_assignment$vax_1[index] == 1) / sum(index)
      
      if (age == "0-17 years") vax_rate <- vax_rate / 0.967
      
      cat(sprintf("%s, %s, %% vaccinated: %.3f\n", age, risk, vax_rate))
    }
  }
  
  # ============================================================================
  # Assign second vaccine dose (vax_2)
  # ============================================================================
  # Only for 65+ and immunocompromised under 65
  # NOTE: Only individuals who have already been assigned with receiving the first vaccine dose `vax_1` are eligible for the second vaccine dose. Since we are sampling only
  #       from individuals who have been assigned first vaccine dose but we are trying to hit the age/risk-specific coverage targets for the entire age/risk group 
  #       (including those who aren't eligible for the second dose), we adjust the coverage upwards for the probability in the rbinom function.
  
  # Define second dose eligible (those who got first dose)
  received_vax_1 <- which(vax_assignment$vax_1 == 1)
  
  assign_vax_2 <- function(this_age_index, this_risk_index, coverage) {
    eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, received_vax_1))
    
    if (length(eligible_index) == 0) return(invisible(NULL))
    
    n_target <- length(Reduce(intersect, list(this_age_index, this_risk_index)))
    
    # Assign vaccination probabilistically
    vax_assignment[["vax_2"]][eligible_index] <<- rbinom(
      n = length(eligible_index), 
      size = 1, 
      prob = coverage * n_target / length(eligible_index)
    )
    
    invisible(NULL)
  }
  
  assign_vax_2(age_18_29_index, risk_immuno_index, 0.054)
  assign_vax_2(age_30_49_index, risk_immuno_index, 0.054)
  assign_vax_2(age_50_64_index, risk_immuno_index, 0.054)
  
  # Age 65-74 years and 75+ (all risk groups)
  assign_vax_2(age_65_74_index, age_65_74_index, 0.089)
  assign_vax_2(age_75_plus_index, age_75_plus_index, 0.089)
  
  # ============================================================================
  # Print statements to check assignments of vax 2
  # ============================================================================
  cat("\n--- Vaccination rates (vax 2) ---\n")
  for (age in c("18-29 years", "30-49 years", "50-64 years")) {
    index <- vax_assignment$age_group == age & vax_assignment$risk_group == "immunocompromised"
    vax_rate <- sum(vax_assignment$vax_2[index] == 1) / sum(index)
    
    cat(sprintf("%s, %s, %% vaccinated: %.3f\n", age, risk, vax_rate))
  }
  
  for (age in c("65-74 years", "75+ years")) {
      index <- vax_assignment$age_group == age 
      vax_rate <- sum(vax_assignment$vax_2[index] == 1) / sum(index)
      
      cat(sprintf("%s, %% vaccinated: %.3f\n", age, vax_rate))
  }
  
  # ============================================================================
  # Assign third vaccine dose (vax_3)
  # ============================================================================
  # Complex logic: some groups get deterministic assignment, others use rbinom
  # Assign who will receive the third vaccine dose (filter first among individuals who have already been assigned with receiving second vaccine dose `vax_2`, then filter among those who received first vaccine dose `vax_1`)
  # NOTE: Only individuals who have been previously vaccinated are eligible for the third vaccine dose. We filter first among individuals who have already been assigned with receiving second vaccine dose `vax_2`, 
  #       then filter among those who received first vaccine dose `vax_1`, then among all previously vaccinated individuals. We adjust the coverage upwards for the probability in the rbinom function.
  
  # received_vax_2 used for those who were eligible for second dose and received it. received_vax_1 used here for those not eligible for second_dose
  # i.e., second dose eligible means they received the 1st dose -- use this for groups not eligible for second dose
  
  received_vax_2 <- which(vax_assignment$vax_2 == 1)
  
  assign_vax_3 <- function(this_age_index, this_risk_index, second_dose_eligible,
                          third_dose_coverage, first_dose_coverage, second_dose_coverage = NULL) {
    
    if (third_dose_coverage < first_dose_coverage) {
      if (second_dose_eligible) {
        
        # Case 1: Assign to everyone who received vax 2 first, then top-up with those who received vax 1
        eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, received_vax_2))
        vax_assignment$vax_3[eligible_index] <<- 1
        
        vax_1_not_2_index <- setdiff(Reduce(intersect, list(this_age_index, this_risk_index, received_vax_1)), received_vax_2)
        n_to_sample <- round(length(Reduce(intersect, list(this_age_index, this_risk_index))) * (third_dose_coverage - second_dose_coverage))
        
        if (n_to_sample > length(vax_1_not_2_index)) {
          print("Error: not enough eligible people to sample from - need to broaden to previously vaccinated")
          break
        }
        
        sampled_index <- sample(vax_1_not_2_index, n_to_sample)
        vax_assignment$vax_3[sampled_index] <<- 1
      } 
      else {
        # Case 2: Sample from those who received vax 1
        eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, received_vax_1))
        n_target <- length(Reduce(intersect, list(this_age_index, this_risk_index)))
        
        vax_assignment$vax_3[eligible_index] <<- rbinom(
          n = length(eligible_index),
          size = 1,
          prob = third_dose_coverage * n_target / length(eligible_index)
        )
      }
    } 
    else if (third_dose_coverage > first_dose_coverage) {
      # Case 3: Assign to everyone who received vax 1, then top-up with those who were previously vaccinated
      eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, received_vax_1))
      vax_assignment$vax_3[eligible_index] <<- 1
      
      prev_not_1_index <- setdiff(
        Reduce(intersect, list(this_age_index, this_risk_index, prev_vaccinated_index)),
        received_vax_1
      )
      
      n_to_sample <- round(length(Reduce(intersect, list(this_age_index, this_risk_index))) * (third_dose_coverage - first_dose_coverage))
        
      if (n_to_sample > length(prev_not_1_index)) {
        print("Error: not enough eligible people to sample from - need to broaden to those not previously vaccinated")
        break
      }
      
      sampled_index <- sample(prev_not_1_index, n_to_sample)
      vax_assignment$vax_3[sampled_index] <<- 1
      }
    
    invisible(NULL)
  }
  
  # Age 0-17 years
  assign_vax_3(age_0_17_index, risk_healthy_index, second_dose_eligible = FALSE, third_dose_coverage = 0.1194, first_dose_coverage = 0.131)
  assign_vax_3(age_0_17_index, risk_highrisk_index, second_dose_eligible = FALSE, third_dose_coverage = 0.2397, first_dose_coverage = 0.235)
  assign_vax_3(age_0_17_index, risk_immuno_index, second_dose_eligible = FALSE, third_dose_coverage = 0.2397, first_dose_coverage = 0.235) 

  # Age 18-29 years
  assign_vax_3(age_18_29_index, risk_healthy_index, second_dose_eligible = FALSE, third_dose_coverage = 0.0758, first_dose_coverage = 0.93)
  assign_vax_3(age_18_29_index, risk_highrisk_index, second_dose_eligible = FALSE, third_dose_coverage = 0.1521, first_dose_coverage = 0.167)
  assign_vax_3(age_18_29_index, risk_immuno_index, second_dose_eligible = TRUE, third_dose_coverage = 0.1521, first_dose_coverage = 0.167, second_dose_coverage = 0.054) 

  # Age 30-49 years
  assign_vax_3(age_30_49_index, risk_healthy_index, second_dose_eligible = FALSE, third_dose_coverage = 0.121, first_dose_coverage = 0.131)
  assign_vax_3(age_30_49_index, risk_highrisk_index, second_dose_eligible = FALSE, third_dose_coverage = 0.2428, first_dose_coverage = 0.2355)
  assign_vax_3(age_30_49_index, risk_immuno_index, second_dose_eligible = TRUE, third_dose_coverage = 0.2428, first_dose_coverage = 0.2355, second_dose_coverage = NULL)
  
  # Age 50-64 years
  assign_vax_3(age_50_64_index, risk_healthy_index, second_dose_eligible = FALSE, third_dose_coverage = 0.1471, first_dose_coverage = 0.1439)
  assign_vax_3(age_50_64_index, risk_highrisk_index, second_dose_eligible = FALSE, third_dose_coverage = 0.2953, first_dose_coverage = 0.2591)
  assign_vax_3(age_50_64_index, risk_immuno_index, second_dose_eligible = TRUE, third_dose_coverage = 0.2953, first_dose_coverage = 0.2591, second_dose_coverage = NULL)
  
  # Age 65-74 years
  assign_vax_3(age_65_74_index, risk_healthy_index, second_dose_eligible = TRUE, third_dose_coverage = 0.2334, first_dose_coverage = 0.2228)
  assign_vax_3(age_65_74_index, risk_highrisk_index, second_dose_eligible = TRUE, third_dose_coverage = 0.4685, first_dose_coverage = 0.4010)
  assign_vax_3(age_65_74_index, risk_immuno_index, second_dose_eligible = TRUE, third_dose_coverage = 0.4685, first_dose_coverage = 0.4010, second_dose_coverage = NULL)
  
  # Age 75+ years
  assign_vax_3(age_75_plus_index, risk_healthy_index, second_dose_eligible = TRUE, third_dose_coverage = 0.254, first_dose_coverage = 0.221)
  assign_vax_3(age_75_plus_index, risk_highrisk_index, second_dose_eligible = TRUE, third_dose_coverage = 0.5097, first_dose_coverage = 0.398)
  assign_vax_3(age_75_plus_index, risk_immuno_index, second_dose_eligible = TRUE, third_dose_coverage = 0.5097, first_dose_coverage = 0.398, second_dose_coverage = NULL)
  
  # ============================================================================
  # Print statements to check assignments of vax 3
  # ============================================================================
  targets <- c("0-17 years" = 0.134, 
               "18-29 years" = 0.097, 
               "30-49 years" = 0.154,
               "50-64 years" = 0.241, 
               "65-74 years" = 0.439, 
               "75+ years" = 0.478)
  
  # Checks for overall age group
  cat("\n--- Vaccination rates by age (vax 3) ---\n")
  
  for (age in age_groups) {
    index <- vax_assignment$age_group == age
    vax_rate <- sum(vax_assignment$vax_3[index] == 1) / sum(index)
    
    # Apply adjustment for 0-17 age group
    if (age == "0-17 years") vax_rate <- vax_rate / 0.967
    
    cat(sprintf("Percentage of %s receiving vaccinated: %.3f (target: %.3f)\n", 
                age, vax_rate, targets[age]))
  }
  
  # Checks for age and risk group
  cat("\n--- Vaccination rates by age and risk group (vax 3) ---\n")
  for (age in age_groups) {
    for (risk in risk_groups) {
      index <- vax_assignment$age_group == age & vax_assignment$risk_group == risk
      vax_rate <- sum(vax_assignment$vax_3[index] == 1) / sum(index)
      
      if (age == "0-17 years") vax_rate <- vax_rate / 0.967
      
      cat(sprintf("%s, %s, %% vaccinated: %.3f\n", age, risk, vax_rate))
    }
  }
  
  # ============================================================================
  # Assign fourth vaccine dose (vax_4)
  # ============================================================================
  # Only for 65+ and immunocompromised under 65
  # NOTE: Sample from those who received vax_2 (all vax_4 uptakes are less than vax_2 -- see Supp. Materials for details)
  
  assign_vax_4 <- function(this_age_index, this_risk_index, coverage) {
    eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, received_vax_2))
    
    if (length(eligible_index) == 0) return(invisible(NULL))
    
    n_target <- length(Reduce(intersect, list(this_age_index, this_risk_index)))
    
    # Assign vaccination probabilistically
    vax_assignment[["vax_4"]][eligible_index] <<- rbinom(
      n = length(eligible_index), 
      size = 1, 
      prob = coverage * n_target / length(eligible_index)
    )
    
    invisible(NULL)
  }
  
  assign_vax_4(age_18_29_index, risk_immuno_index, 0.033)
  assign_vax_4(age_30_49_index, risk_immuno_index, 0.033)
  assign_vax_4(age_50_64_index, risk_immuno_index, 0.033)
  
  # Age 65-74 years and 75+ (all risk groups)
  assign_vax_4(age_65_74_index, age_65_74_index, 0.057)
  assign_vax_4(age_75_plus_index, age_75_plus_index, 0.057)
  
  # ============================================================================
  # Print statements to check assignments of vax 4
  # ============================================================================
  cat("\n--- Vaccination rates (vax 4) ---\n")
  for (age in c("18-29 years", "30-49 years", "50-64 years")) {
    index <- vax_assignment$age_group == age & vax_assignment$risk_group == "immunocompromised"
    vax_rate <- sum(vax_assignment$vax_4[index] == 1) / sum(index)

    cat(sprintf("%s, %s, %% vaccinated: %.3f\n", age, risk, vax_rate))
  }

  for (age in c("65-74 years", "75+ years")) {
    index <- vax_assignment$age_group == age
    vax_rate <- sum(vax_assignment$vax_4[index] == 1) / sum(index)

    cat(sprintf("%s, %% vaccinated: %.3f\n", age, vax_rate))
  }
  
  # ============================================================================
  # Assign vaccine timing - First dose (days 1-182)
  #     Distribute first vaccine dose over first 365 days, according to age-specific distribution of historic coverage
  # ============================================================================
  
  assign_vax_timing <- function(timing_col,          # vaccine_wave1, vaccine_wave2, vaccine_wave3
                                dist_time,           # time period during which vaccines are distributed (e.g., c(1:365))
                                to_vaccinate_index,  # e.g., which(vax_assignment$vax_1 == 1)
                                timing_data,         # e.g., dose1_data
                                this_age_group,      # e.g., "18-29 years"
                                this_risk_group      # e.g., "immunocompromised"
                                ) {
    
    this_age_index <- age_mapping[[this_age_group]]
    this_risk_index <- risk_mapping[[this_risk_group]]
    
    # Get indices for this age-risk-vaccination status combination
    eligible_index <- Reduce(intersect, list(this_age_index, this_risk_index, to_vaccinate_index))
    
    # Skip if no one to vaccinate
    if (length(eligible_index) == 0) return(invisible(NULL))
    
    # Get timing probabilities from data
    timing_probs <- timing_data %>%
      filter(age_group == this_age_group, risk_group == this_risk_group) %>%
      arrange(week) %>%
      pull(prop_up_to_date_sim_period)
    
    # Create daily probabilities: first day gets first value, then repeat each week value 7 times
    if (timing_col == "vaccine_wave1") {
      daily_probs <- c(timing_probs[1], rep(timing_probs, each = 7))
    } else if (timing_col == "vaccine_wave4") {
      daily_probs <- c(timing_probs[1], rep(timing_probs[2:length(timing_probs)], each = 7))
    } else {
      daily_probs <- c(rep(timing_probs, each = 7))
    }
    
    # Assign vaccination days
    vax_assignment[[timing_col]][eligible_index] <<- sample(
      x = dist_time,
      size = length(eligible_index),
      prob = daily_probs,
      replace = TRUE
    )
    
    invisible(NULL)
  }
  
  age_mapping <- list(
    "0-17 years" = age_0_17_index,
    "18-29 years" = age_18_29_index,
    "30-49 years" = age_30_49_index,
    "50-64 years" = age_50_64_index,
    "65-74 years" = age_65_74_index,
    "75+ years" = age_75_plus_index
  )
  
  risk_mapping <- list(
    "healthy" = risk_healthy_index,
    "higher risk" = risk_highrisk_index,
    "immunocompromised" = risk_immuno_index
  )
  
  # Assign first dose timing for all age-risk combinations
  for (age_name in names(age_mapping)) {
    for (risk_name in names(risk_mapping)) {
      assign_vax_timing(
        timing_col = "vaccine_wave1",
        dist_time = c(1:365),
        to_vaccinate_index = which(vax_assignment$vax_1 == 1),
        timing_data = dose1_data,
        this_age_group = age_name,
        this_risk_group = risk_name
      )
    }
  }
  
  # Assign second dose timing for all age-risk combinations
  for (age_name in names(age_mapping)) {
    for (risk_name in names(risk_mapping)) {
      assign_vax_timing(
        timing_col = "vaccine_wave2",
        dist_time = c(240:365),
        to_vaccinate_index = which(vax_assignment$vax_2 == 1),
        timing_data = dose2_data,
        this_age_group = age_name,
        this_risk_group = risk_name
      )
    }
  }
  
  # Assign third dose timing for all age-risk combinations
  for (age_name in names(age_mapping)) {
    for (risk_name in names(risk_mapping)) {
      assign_vax_timing(
        timing_col = "vaccine_wave3",
        dist_time = c(366:547),
        to_vaccinate_index = which(vax_assignment$vax_3 == 1),
        timing_data = dose3_data,
        this_age_group = age_name,
        this_risk_group = risk_name
      )
    }
  }
  
  # Assign fourth dose timing for all age-risk combinations
  for (age_name in names(age_mapping)) {
    for (risk_name in names(risk_mapping)) {
      assign_vax_timing(
        timing_col = "vaccine_wave4",
        dist_time = c(609:665),
        to_vaccinate_index = which(vax_assignment$vax_4 == 1),
        timing_data = dose4_data,
        this_age_group = age_name,
        this_risk_group = risk_name
      )
    }
  }
  
  return(vax_assignment %>% arrange(individual))
}


