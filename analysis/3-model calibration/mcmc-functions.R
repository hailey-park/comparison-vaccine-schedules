########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

prediction <- function (params) {

  #Read parameters needed for model setup
  baseline_case_hosp_frac <<- params[19]
  time_since <<- params[20]
  
  #Set-up model inputs with updated parameters
  source(here::here("model-setup-eachtime-1mil.R"))
  
  #Run model
  set.seed(88)
  results <- historicalVaccinationSimulation(clean_df[[1]], c(params, lambdas)) 
  
  #Reformat results into weekly incidence estimates by age group (separate columns)
  reformatted_results <- melt(setDT(results %>% group_by(age_group) %>% summarise(across(total_pop:day547, sum))), id = c("age_group", "total_pop"))  %>%
    mutate(days = as.numeric(str_sub(variable, 4)),
           weeks = ceiling(days/7)) %>%
    group_by(age_group, weeks) %>% 
    summarise(total_cases = sum(value), total_pop = mean(total_pop)) %>%
    mutate(simulated_inc = total_cases/total_pop * 100000) %>%
    dplyr::select(-c(total_pop, total_cases)) %>%
    group_by(age_group) %>%
    pivot_wider(names_from = age_group, values_from = simulated_inc) %>%
    dplyr::select(-`0-17 years`)  %>% #removing 0-17 years (severe incidence is set to 0 so no need to fit)
    filter(weeks != 79) #remove 79rd week of results (only captures 365th day)

    return(reformatted_results)
}


# The score is the dynamic time warping (DTW) distance between the simulated inc predictions 
# and the observed data (units = weeks).
score <- function (sim_pred, data) {

  age_18_29 <- dtw(data$`18-29 years`, sim_pred$`18-29 years`,
                  window.type = "sakoechiba",
                  window.size = 2,
                  keep=TRUE)$distance
  age_30_49 <- dtw(data$`30-49 years`, sim_pred$`30-49 years`,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_50_64 <- dtw(data$`50-64 years`, sim_pred$`50-64 years`,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_65_74 <- dtw(data$`65-74 years`, sim_pred$`65-74 years`,
                   window.type = "sakoechiba",
                   window.size = 2,
                   keep=TRUE)$distance
  age_75_plus <- dtw(data$`75+ years`, sim_pred$`75+ years`,
                     window.type = "sakoechiba",
                     window.size = 2,
                     keep=TRUE)$distance

  total <- sum(age_18_29, age_30_49, age_50_64, age_65_74, age_75_plus)
  return(total)
}


#This is a stricter (narrower) sigmoid function used for the acceptance ratio
sigmoid <- function(x, midpoint = -30, max_y = 10000, min_y = 0, max_diff_x = 50, max_diff_y = 300){
  k <- ((max_y - max_diff_y)/(max_diff_y - min_y))^(1/(max_diff_x - midpoint))
  y <- (max_y * k^midpoint + min_y * k^x) / (k^midpoint  + k^x)
  y_rescaled <- y/10000
  return(y_rescaled)
}

# Typically in MCMC, the likelihood is the probability (density) with which we would expect the 
# observed data to occur conditional on the parameters of the model that we look at.
# For this MCMC, we just have the likelihood as DTW scoring.
likelihood <- function (params) {
  
  #Get predictions
  sim_pred <- prediction(params)
  return(score(sim_pred, observed_data))
}


#Typically in MCMC, the posterior distribution would be some combination of the likelihood of the proposed parameter
#set with the distribution of priors, but since our parameter seach is uninformed by priors, posterior distribution
#is just likelihood.
posterior <- function(param){
  likel <- likelihood(param)
  print(paste0("Likelihood: ", likel))
  return (likel)
  }


# Choosing a new parameter value close to the old value based on some 
# probability density that is called the proposal function. Here, the 
# proposal is a constrained normal distribution centered at the current value
proposalfunction <- function(param){
  
  lambda_1 = min(max(rnorm(1, mean=param[1], sd=0.0005), 0.98), 1.01)
  lambda_2 = min(max(rnorm(1, mean=param[2], sd=0.0005), 0.98), 1.01)
  lambda_3 = min(max(rnorm(1, mean=param[3], sd=0.0005), 0.98), 1.01)
  lambda_4 = min(max(rnorm(1, mean=param[4], sd=0.0005), 0.98), 1.01)
  lambda_5 = min(max(rnorm(1, mean=param[5], sd=0.0005), 0.98), 1.01)
  lambda_6 = min(max(rnorm(1, mean=param[6], sd=0.0005), 0.98), 1.01)
  lambda_7 = min(max(rnorm(1, mean=param[7], sd=0.0005), 0.98), 1.01)
  lambda_8 = min(max(rnorm(1, mean=param[8], sd=0.0005), 0.98), 1.01)
  lambda_9 = min(max(rnorm(1, mean=param[9], sd=0.0005), 0.98), 1.01)
  lambda_10 = min(max(rnorm(1, mean=param[10], sd=0.0005), 0.98), 1.01)
  lambda_11 = min(max(rnorm(1, mean=param[11], sd=0.0005), 0.98), 1.01)
  lambda_12 = min(max(rnorm(1, mean=param[12], sd=0.0005), 0.98), 1.01)
  lambda_13 = min(max(rnorm(1, mean=param[13], sd=0.0005), 0.98), 1.01)
  lambda_14 = min(max(rnorm(1, mean=param[14], sd=0.0005), 0.98), 1.01)
  lambda_15 = min(max(rnorm(1, mean=param[15], sd=0.0005), 0.98), 1.01)
  lambda_16 = min(max(rnorm(1, mean=param[16], sd=0.0005), 0.98), 1.01)
  lambda_17 = min(max(rnorm(1, mean=param[17], sd=0.0005), 0.98), 1.01)
  lambda_18 = min(max(rnorm(1, mean=param[18], sd=0.0005), 0.98), 1.01)
  
  baseline_case_hosp_frac = min(max(rnorm(1, mean=param[19], sd=param[19] * 0.1), 0.001), 0.05)
  time_since = min(max(rnorm(1, mean=param[20], sd=1), -10), 10)
  
  return(c(lambda_1, lambda_2, lambda_3, lambda_4, lambda_5, lambda_6, lambda_7, 
           lambda_8, lambda_9, lambda_10, lambda_11, lambda_12, lambda_13,
           lambda_14, lambda_15, lambda_16, lambda_17, lambda_18,
           baseline_case_hosp_frac, time_since))
}


run_metropolis_MCMC <- function(startvalue, max_iteration, starting_index){
  
  #Either read in existing chain, or start new one
  if (starting_index != 1) {
    chain = as.matrix(read.csv("mcmc-output/mcmc-params-18lambdas-window2-052325.csv")[,-1])
  } else{
    chain = array(dim = c(max_iteration+1, 20 + 2))
    chain[1,] = c(startvalue, posterior(startvalue), NA)
  }
  
  #Run mcmc for `max_iteration`
  for (i in starting_index:max_iteration){
    print(i)
    print(chain[i,1:20])
    
    # Generate a new candidate sample from the proposal distribution
    proposal = proposalfunction(chain[i,1:20])
    print("Proposal: ")
    print(proposal)
    
    # Compute posterior probability of proposal vs current state
    posterior_proposal <- posterior(proposal)
    posterior_current <- chain[i,21]
    
    # Compute acceptance ratio (using sigmoid scaled to (0,1))
    # If random draw is less than acceptance ratio, accept the proposal
    acceptance_ratio <- sigmoid((posterior_proposal - posterior_current))
    random_probab <- runif(1)
    
    print(paste0("Probability: ", acceptance_ratio))
    print(paste0("Random prob: ", random_probab))
    
    if (random_probab < acceptance_ratio){
      chain[i+1,1:20] <- proposal
      chain[i+1,21] <- posterior_proposal
      chain[i+1,22] <- acceptance_ratio
      
    } else{
      chain[i+1,1:20] <- chain[i, 1:20]
      chain[i+1,21] <- posterior_current
      chain[i+1,22] <- acceptance_ratio
      
    }
    
    print(paste0("Iteration ", i, " Chain: "))
    print(chain[i+1,])
    
    # Every simulation, write the chain to .csv (since each simulation run takes a few min so calibration takes much longer than couple of hours)
    if (i%%1 == 0) {
      write.csv(chain, "mcmc-output/mcmc-params-18lambdas-window2-052325.csv")
    }
  }
  return(chain)
}
