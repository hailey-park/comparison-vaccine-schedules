########################################################################################################################
# Title: Vaccination assignment under historical coverage
# Author: Hailey Park
# Date: January 16th, 2025
# Last updated: January 28, 2026
########################################################################################################################
# OVERVIEW:
# This function assigns 
#       - vaccine doses (vax_1, vax_2, vax_3) 
#       - timing (vaccine_wave, second_vaccine_wave, third_vaccine_wave) 
# based on historical age/risk-specific vaccine administration data.
# 
# Semi-annual doses (vax_2) are only assigned to 65+ years and 
# immunocompromised groups based on coverage estimates.
#
# INPUTS:
#   df: data frame number of people x 11 columns that includes prior immune status, age, and risk group
#   first, second, third, and fourth, dose data: data frames with corresponding vaccine coverage by age and risk over time
#
# OUTPUT:
#   vax_assignment: same as input df with additional columns for vax_1, vax_2, vax_3 and timing
#

historical_vax_assignment <- function(df, first_dose_data, second_dose_data, next_year_dose_data){
    
  # Initialize vaccination assignment columns
    vax_assignment <- df %>% 
      mutate(vax_1 = 0, vax_2 = 0, vax_3 = 0, vaccine_wave = 0, second_vaccine_wave = 0, third_vaccine_wave = 0)
    
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
    prev_vaccinated <- which(vax_assignment$prior_vacc != "unvax")
    
    # ============================================================================
    # Assign first vaccine dose (vax_1)
    # ============================================================================
    # Uses age-specific and risk-specific coverage rates
    # Only previously vaccinated individuals are eligible
    
    # Age 0-17 years
    vax_assignment$vax_1[Reduce(intersect, list(age_0_17_index,risk_healthy_index,prev_vaccinated))] <- 
      rbinom(length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,prev_vaccinated))), 1, 
             0.13 * length(Reduce(intersect, list(age_0_17_index,risk_healthy_index))) / 
               length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,prev_vaccinated))))
    
    vax_assignment$vax_1[Reduce(intersect, list(age_0_17_index,risk_highrisk_index,prev_vaccinated))] <- 
      rbinom(length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index,prev_vaccinated))), 1, 
             0.235 * length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index))) / 
               length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index,prev_vaccinated))))
    
    vax_assignment$vax_1[Reduce(intersect, list(age_0_17_index,risk_immuno_index,prev_vaccinated))] <- 
      rbinom(length(Reduce(intersect, list(age_0_17_index,risk_immuno_index,prev_vaccinated))), 1, 
             0.235 * length(Reduce(intersect, list(age_0_17_index,risk_immuno_index))) / 
               length(Reduce(intersect, list(age_0_17_index,risk_immuno_index,prev_vaccinated))))
    
    # Age 18-29 years
    vax_assignment$vax_1[Reduce(intersect, list(age_18_29_index,risk_healthy_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,prev_vaccinated))), 1, 0.0926 * length(Reduce(intersect, list(age_18_29_index,risk_healthy_index)))/length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_18_29_index,risk_highrisk_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,prev_vaccinated))), 1, 0.1668 * length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_18_29_index,risk_immuno_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,prev_vaccinated))), 1, 0.1668 * length(Reduce(intersect, list(age_18_29_index,risk_immuno_index)))/length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,prev_vaccinated))))
    
    # Age 30-49 years
    vax_assignment$vax_1[Reduce(intersect, list(age_30_49_index,risk_healthy_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,prev_vaccinated))), 1, 0.1308 * length(Reduce(intersect, list(age_30_49_index,risk_healthy_index)))/length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_30_49_index,risk_highrisk_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index,prev_vaccinated))), 1, 0.2355 * length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_30_49_index,risk_immuno_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,prev_vaccinated))), 1, 0.2355 * length(Reduce(intersect, list(age_30_49_index,risk_immuno_index)))/length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,prev_vaccinated))))
    
    # Age 50-64 years
    vax_assignment$vax_1[Reduce(intersect, list(age_50_64_index,risk_healthy_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_50_64_index,risk_healthy_index,prev_vaccinated))), 1, 0.1439 * length(Reduce(intersect, list(age_50_64_index,risk_healthy_index)))/length(Reduce(intersect, list(age_50_64_index,risk_healthy_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_50_64_index,risk_highrisk_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index,prev_vaccinated))), 1, 0.2591 * length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_50_64_index,risk_immuno_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,prev_vaccinated))), 1, 0.2591 * length(Reduce(intersect, list(age_50_64_index,risk_immuno_index)))/length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,prev_vaccinated))))
    
    # Age 65-74 years
    vax_assignment$vax_1[Reduce(intersect, list(age_65_74_index,risk_healthy_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_65_74_index,risk_healthy_index,prev_vaccinated))), 1, 0.2228 * length(Reduce(intersect, list(age_65_74_index,risk_healthy_index)))/length(Reduce(intersect, list(age_65_74_index,risk_healthy_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_65_74_index,risk_highrisk_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index,prev_vaccinated))), 1, 0.4010 * length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_65_74_index,risk_immuno_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_65_74_index,risk_immuno_index,prev_vaccinated))), 1, 0.4010 * length(Reduce(intersect, list(age_65_74_index,risk_immuno_index)))/length(Reduce(intersect, list(age_65_74_index,risk_immuno_index,prev_vaccinated))))
    
    # Age 75+ years
    vax_assignment$vax_1[Reduce(intersect, list(age_75_plus_index,risk_healthy_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index,prev_vaccinated))), 1, 0.221 * length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index)))/length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,prev_vaccinated))), 1, 0.3979 * length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,prev_vaccinated))))
    vax_assignment$vax_1[Reduce(intersect, list(age_75_plus_index,risk_immuno_index,prev_vaccinated))] <- rbinom(length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index,prev_vaccinated))), 1, 0.3979 * length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index)))/length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index,prev_vaccinated))))
    
    print(paste0("Percentage of 0-17 years receiving vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$vax_1 == 1))/(length(which(vax_assignment$age_group == "0-17 years"))* 0.967)))
    print(paste0("Percentage of 18-29 years receiving vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "18-29 years"))))
    print(paste0("Percentage of 30-49 years receiving vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "30-49 years"))))
    print(paste0("Percentage of 50-64 years receiving vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "50-64 years"))))
    print(paste0("Percentage of 65-74 years receiving vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "65-74 years"))))
    print(paste0("Percentage of 75+ years receiving vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "75+ years"))))
    print("0.15, 0.113, 0.1596, 0.217, 0.3794, 0.3764")
    
    print(paste0("0-17 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "healthy"))* 0.967)))
    print(paste0("18-29 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("30-49 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("50-64 years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("65-74 years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("75+ years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "healthy"))))
    
    print(paste0("0-17 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "higher risk"))* 0.967)))
    print(paste0("18-29 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("30-49 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("50-64 years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("65-74 years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("75+ years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "higher risk"))))
    
    print(paste0("0-17 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "immunocompromised"))* 0.967)))
    print(paste0("18-29 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("30-49 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("50-64 years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("65-74 years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("75+ years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_1 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "immunocompromised"))))
    
    # ============================================================================
    # Assign second vaccine dose (vax_2)
    # ============================================================================
    # Only for 65+ and immunocompromised under 65
    # NOTE: Only individuals who have already been assigned with receiving the first vaccine dose `vax_1` are eligible for the second vaccine dose. Since we are sampling only
    #       from individuals who have been assigned first vaccine dose but we are trying to hit the age/risk-specific coverage targets for the entire age/risk group 
    #       (including those who aren't eligible for the second dose), we adjust the coverage upwards for the probability in the rbinom function.
    
    # Define second dose eligible (those who got first dose)
    second_dose_eligible <- which(vax_assignment$vax_1 == 1)
    
    # Immunocompromised under 65
    vax_assignment$vax_2[Reduce(intersect, list(age_18_29_index,risk_immuno_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,second_dose_eligible))), 1, 0.054 * length(Reduce(intersect, list(age_18_29_index,risk_immuno_index)))/length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,second_dose_eligible))))
    vax_assignment$vax_2[Reduce(intersect, list(age_30_49_index,risk_immuno_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,second_dose_eligible))), 1, 0.054 * length(Reduce(intersect, list(age_30_49_index,risk_immuno_index)))/length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,second_dose_eligible))))
    vax_assignment$vax_2[Reduce(intersect, list(age_50_64_index,risk_immuno_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,second_dose_eligible))), 1, 0.054 * length(Reduce(intersect, list(age_50_64_index,risk_immuno_index)))/length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,second_dose_eligible))))
    
    # Age 65-74 years and 75+ (all risk groups)
    vax_assignment$vax_2[Reduce(intersect, list(age_65_74_index, second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_65_74_index, second_dose_eligible))), 1, 0.089 * length(Reduce(intersect, list(age_65_74_index)))/length(Reduce(intersect, list(age_65_74_index,second_dose_eligible))))
    vax_assignment$vax_2[Reduce(intersect, list(age_75_plus_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_75_plus_index, second_dose_eligible))), 1, 0.089 * length(Reduce(intersect, list(age_75_plus_index)))/length(Reduce(intersect, list(age_75_plus_index,second_dose_eligible))))
    
    print(paste0("Percentage of immunocompromised receiving second dose: ",  length(which(vax_assignment$risk_group %in% c("immunocompromised") & vax_assignment$vax_2 == 1))/length(which(vax_assignment$risk_group %in% c("immunocompromised")))))
    print(paste0("Percentage of 65-74 years receiving second dose: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$vax_2 == 1))/length(which(vax_assignment$age_group == "65-74 years"))))
    print(paste0("Percentage of 75+ years receiving second dose: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$vax_2 == 1))/length(which(vax_assignment$age_group == "75+ years"))))
    
    # ============================================================================
    # Assign third vaccine dose (vax_3)
    # ============================================================================
    # Complex logic: some groups get deterministic assignment, others use rbinom
    #Assign who will receive the third vaccine dose (filter first among individuals who have already been assigned with receiving second vaccine dose `vax_2`, then filter among those who received first vaccine dose `vax_1`)
    # NOTE: Only individuals who have been previously vaccinated are eligible for the third vaccine dose. We filter first among individuals who have already been assigned with receiving second vaccine dose `vax_2`, 
    #       then filter among those who received first vaccine dose `vax_1`, then among all previously vaccinated individuals. We adjust the coverage upwards for the probability in the rbinom function.
    third_dose_eligible <- which(vax_assignment$vax_2 == 1)
    
    # KMB: should the term here be third_dose_eligible instead of second_dose_eligible?
    # second dose eligible means they received the 1st dose -- use this for groups not eligible for second dose
    
    # Age 0-17 years
    vax_assignment$vax_3[Reduce(intersect, list(age_0_17_index,risk_healthy_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,second_dose_eligible))), 1, 0.1194 * length(Reduce(intersect, list(age_0_17_index,risk_healthy_index)))/length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,second_dose_eligible))))
    vax_assignment$vax_3[Reduce(intersect, list(age_0_17_index,risk_highrisk_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_0_17_index,risk_highrisk_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index))) * (0.2397 - 0.235)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_0_17_index,risk_immuno_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_0_17_index,risk_immuno_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_0_17_index,risk_immuno_index))) *  (0.2397 - 0.235)))] <- 1
    
    # Age 18-29 years
    vax_assignment$vax_3[Reduce(intersect, list(age_18_29_index,risk_healthy_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,second_dose_eligible))), 1, 0.0758 * length(Reduce(intersect, list(age_18_29_index,risk_healthy_index)))/length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,second_dose_eligible))))
    vax_assignment$vax_3[Reduce(intersect, list(age_18_29_index,risk_highrisk_index,second_dose_eligible))] <-  rbinom(length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,second_dose_eligible))), 1, 0.1521 * length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index)))/length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,second_dose_eligible))))
    vax_assignment$vax_3[Reduce(intersect, list(age_18_29_index,risk_immuno_index,third_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_18_29_index,risk_immuno_index, second_dose_eligible)), third_dose_eligible), round(length(Reduce(intersect, list(age_18_29_index,risk_immuno_index))) *  (0.1521 - 0.054)))] <- 1
    
    # Age 30-49 years
    vax_assignment$vax_3[Reduce(intersect, list(age_30_49_index,risk_healthy_index,second_dose_eligible))] <- rbinom(length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,second_dose_eligible))), 1, 0.121 * length(Reduce(intersect, list(age_30_49_index,risk_healthy_index)))/length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,second_dose_eligible))))
    vax_assignment$vax_3[Reduce(intersect, list(age_30_49_index,risk_highrisk_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_30_49_index,risk_highrisk_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index))) * (0.2428 - 0.2355)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_30_49_index,risk_immuno_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_30_49_index,risk_immuno_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_30_49_index,risk_immuno_index))) * (0.2428 - 0.2355)))] <- 1
    
    # Age 50-64 years
    vax_assignment$vax_3[Reduce(intersect, list(age_50_64_index,risk_healthy_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_50_64_index,risk_healthy_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_50_64_index,risk_healthy_index))) * (0.1471 - 0.1439)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_50_64_index,risk_highrisk_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_50_64_index,risk_highrisk_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index))) * (0.2953 - 0.2591)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_50_64_index,risk_immuno_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_50_64_index,risk_immuno_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_50_64_index,risk_immuno_index))) * (0.2953 - 0.2591)))] <- 1
    
    # Age 65-74 years
    vax_assignment$vax_3[Reduce(intersect, list(age_65_74_index,risk_healthy_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_65_74_index,risk_healthy_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_65_74_index,risk_healthy_index))) * (0.2334 - 0.2228)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_65_74_index,risk_highrisk_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_65_74_index,risk_highrisk_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index))) * (0.4685 - 0.4010)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_65_74_index,risk_immuno_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_65_74_index,risk_immuno_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_65_74_index,risk_immuno_index))) * (0.4685 - 0.4010)))] <- 1
    
    # Age 75+ years
    vax_assignment$vax_3[Reduce(intersect, list(age_75_plus_index,risk_healthy_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_75_plus_index,risk_healthy_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index))) * (0.254 - 0.221)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index))) * (0.5097 - 0.398)))] <- 1
    vax_assignment$vax_3[Reduce(intersect, list(age_75_plus_index,risk_immuno_index,second_dose_eligible))] <- 1
    vax_assignment$vax_3[sample(setdiff(Reduce(intersect, list(age_75_plus_index,risk_immuno_index, prev_vaccinated)), second_dose_eligible), round(length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index))) * (0.5097 - 0.398)))] <- 1
    
    # Print coverage summary
    print(paste0("Percentage of 0-17 years receiving vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$vax_3 == 1))/(length(which(vax_assignment$age_group == "0-17 years"))* 0.967)))
    print(paste0("Percentage of 18-29 years receiving vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "18-29 years"))))
    print(paste0("Percentage of 30-49 years receiving vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "30-49 years"))))
    print(paste0("Percentage of 50-64 years receiving vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "50-64 years"))))
    print(paste0("Percentage of 65-74 years receiving vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "65-74 years"))))
    print(paste0("Percentage of 75+ years receiving vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "75+ years"))))
    
    print(paste0("0-17 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "healthy")) * 0.967)))
    print(paste0("18-29 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("30-49 years, healthy: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("50-64 years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("65-74 years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "healthy"))))
    print(paste0("75+ years, healthy: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "healthy" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "healthy"))))
    print("0.1336, 0.0967, 0.1543, 0.2411, 0.4393, 0.478")
    
    print(paste0("0-17 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "higher risk"))* 0.967)))
    print(paste0("18-29 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("30-49 years, higher risk: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("50-64 years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("65-74 years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "higher risk"))))
    print(paste0("75+ years, higher risk: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "higher risk" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "higher risk"))))
    
    print(paste0("0-17 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/(length(which(vax_assignment$age_group == "0-17 years" & vax_assignment$risk_group == "immunocompromised"))* 0.967)))
    print(paste0("18-29 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "18-29 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("30-49 years, immunocompromised: % vaccinated: ", length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "30-49 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("50-64 years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "50-64 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("65-74 years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "65-74 years" & vax_assignment$risk_group == "immunocompromised"))))
    print(paste0("75+ years, immunocompromised: % vaccinated: ",  length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "immunocompromised" & vax_assignment$vax_3 == 1))/length(which(vax_assignment$age_group == "75+ years" & vax_assignment$risk_group == "immunocompromised"))))
    
    # ============================================================================
    # Assign vaccine timing - First dose (days 1-182)
    # ============================================================================
    #distribute first vaccine dose over first 365 days, according to age-specific distribution of historic coverage
    
    to_vaccinate_index <- which(vax_assignment$vax_1 == 1)
    
    # Age 0-17 years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_0_17_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "0-17 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 18-29 years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "18-29 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 30-49 years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "30-49 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 50-64 years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_50_64_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "50-64 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 65-74 years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "65-74 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 75+ years - first dose timing
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "75+ years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "75+ years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "75+ years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "75+ years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_index))] <- sample(c(1:365), length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_index))), prob = c(((first_dose_data %>% filter(age_group == "75+ years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period)[1], rep((first_dose_data %>% filter(age_group == "75+ years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # ============================================================================
    # Assign vaccine timing - Second dose (days 240-365)
    # ============================================================================
    #distribute second vaccine dose over 240-365 days of simulation, according to age-specific distribution of historic coverage
    
    to_vaccinate_second_index <- which(vax_assignment$vax_2 == 1)
    
    # Immunocompromised under 65 - second dose timing
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "18-29 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "30-49 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "50-64 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 65-74 years - second dose timing
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "65-74 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "65-74 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "65-74 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 75+ years - second dose timing
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "75+ years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "75+ years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$second_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_second_index))] <- sample(c(240:365), length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_second_index))), prob = c(rep((second_dose_data %>% filter(age_group == "75+ years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # ============================================================================
    # Assign vaccine timing - Third dose (days 366-547)
    # ============================================================================
    #distribute third vaccine dose last 182 days of simulation, according to age-specific distribution of historic coverage
    to_vaccinate_third_index <- which(vax_assignment$vax_3 == 1)
    
    # Age 0-17 years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_0_17_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "0-17 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_0_17_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "0-17 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_0_17_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_0_17_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "0-17 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 18-29 years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_18_29_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "18-29 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_18_29_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "18-29 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_18_29_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "18-29 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 30-49 years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_30_49_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "30-49 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_30_49_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "30-49 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_30_49_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "30-49 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 50-64 years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_50_64_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "50-64 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_50_64_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "50-64 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_50_64_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "50-64 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 65-74 years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_65_74_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "65-74 years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_65_74_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "65-74 years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_65_74_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "65-74 years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    # Age 75+ years - third dose timing
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_75_plus_index,risk_healthy_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "75+ years", risk_group == "healthy") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_75_plus_index,risk_highrisk_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "75+ years", risk_group == "higher risk") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    vax_assignment$third_vaccine_wave[Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_third_index))] <- sample(c((1 + 365):(182 + 365)), length(Reduce(intersect, list(age_75_plus_index,risk_immuno_index,to_vaccinate_third_index))), prob = c(rep((next_year_dose_data %>% filter(age_group == "75+ years", risk_group == "immunocompromised") %>% arrange(week))$prop_up_to_date_sim_period, each = 7)), replace = TRUE)
    
    return(vax_assignment %>% arrange(individual))
}


