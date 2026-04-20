# Functions of optimisation of natural history parameters - 4-stage version

library(expm)  # for matrix exponentiation

# Function to create transition matrix for 4-stage model
# states:
# 1: undiagnosed stage 1
# 2: undiagnosed stage 2
# 3: undiagnosed stage 3
# 4: undiagnosed stage 4
# 5: diagnosed stage 1
# 6: diagnosed stage 2
# 7: diagnosed stage 3
# 8: diagnosed stage 4

create_transition_matrix <- function(tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                     multiplier = 0 # 1 for multiplicative
                                     ) {
  # tp12: undiagnosed stage 1 -> undiagnosed stage 2
  # tp23: undiagnosed stage 2 -> undiagnosed stage 3
  # tp34: undiagnosed stage 3 -> undiagnosed stage 4
  # tp1: undiagnosed stage 1 -> diagnosed stage 1
  # tp2: undiagnosed stage 2 -> diagnosed stage 2
  # tp3: undiagnosed stage 3 -> diagnosed stage 3
  # tp4: undiagnosed stage 4 -> diagnosed stage 4
  
  ###
  # add effects of age, sex and ethnicity on diagnostic TPs
  
  # tp1 <- tp1*multiplier
  # tp2 <- tp2*multiplier
  # tp3 <- tp3*multiplier
  # tp4 <- tp4*multiplier
  
  # update: add effect on hazard
  # Assume tp4 is unaffected by age, sex and ethnicity
  # lambda_adj <- -log(1 - c(tp1, tp2, tp3)) * multiplier
  # tp_adj <- 1 - exp(-lambda_adj)
  # tp1 <- tp_adj[1]
  # tp2 <- tp_adj[2]
  # tp3 <- tp_adj[3] 
  
  # update: use logit 
  logit_tp <- log(c(tp1, tp2, tp3) / (1 - c(tp1, tp2, tp3)))
  logit_tp_adj <- logit_tp + multiplier
  tp_adj <- exp(logit_tp_adj) / (1 + exp(logit_tp_adj))
  tp1 <- tp_adj[1]
  tp2 <- tp_adj[2]
  tp3 <- tp_adj[3]
  ###
  
  # Cap diagnosis probabilities at maximum feasible values
  # Maximum tp1 = 1 - tp12 -0.001 (to ensure stay_stage1 > 0), and so on
  tp1 <- min(tp1, 1 - tp12 - 0.001)
  tp2 <- min(tp2, 1 - tp23 - 0.001)
  tp3 <- min(tp3, 1 - tp34 - 0.001)
  tp4 <- min(tp4, 1)
  
  # Probability of staying in each undiagnosed stage
  stay_stage1 <- 1 - tp12 - tp1
  stay_stage2 <- 1 - tp23 - tp2
  stay_stage3 <- 1 - tp34 - tp3
  stay_stage4 <- 1 - tp4
  
  # Create 8x8 transition matrix
  P <- matrix(0, nrow = 8, ncol = 8)
  
  # From undiagnosed stage 1 (state 1)
  P[1, 1] <- stay_stage1  # stay in undiagnosed stage 1
  P[1, 2] <- tp12         # progress to undiagnosed stage 2
  P[1, 5] <- tp1          # diagnosed at stage 1
  
  # From undiagnosed stage 2 (state 2)
  P[2, 2] <- stay_stage2  # stay in undiagnosed stage 2
  P[2, 3] <- tp23         # progress to undiagnosed stage 3
  P[2, 6] <- tp2          # diagnosed at stage 2
  
  # From undiagnosed stage 3 (state 3)
  P[3, 3] <- stay_stage3  # stay in undiagnosed stage 3
  P[3, 4] <- tp34         # progress to undiagnosed stage 4
  P[3, 7] <- tp3          # diagnosed at stage 3
  
  # From undiagnosed stage 4 (state 4)
  P[4, 4] <- stay_stage4  # stay in undiagnosed stage 4
  P[4, 8] <- tp4          # diagnosed at stage 4
  
  # Diagnosed states are absorbing (states 5-8)
  P[5, 5] <- 1  # diagnosed stage 1
  P[6, 6] <- 1  # diagnosed stage 2
  P[7, 7] <- 1  # diagnosed stage 3
  P[8, 8] <- 1  # diagnosed stage 4
  
  return(P)
}


# Function to simulate Markov chain for given initial state
simulate_markov <- function(P, initial_state, n_cycles = 24) {
  # Initialize state distribution
  state_dist <- rep(0, 8)
  state_dist[initial_state] <- 1
  
  # Store results for each cycle
  results <- matrix(0, nrow = n_cycles + 1, ncol = 8)
  results[1, ] <- state_dist
  
  # Simulate each cycle
  for (cycle in 1:n_cycles) {
    state_dist <- state_dist %*% P
    results[cycle + 1, ] <- state_dist
  }
  
  return(results)
}


# Function to calculate stage distribution at entry using simulation
# This replaces the analytical approach from the 2-stage model
initial_stage_distribution <- function(month_diagnosis, diagnosed_stage, 
                                       tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                       multiplier = 0, # 1,
                                       n_sim = 10000) {

  ###
  # add effects of age, sex and ethnicity on diagnostic TPs
  
  # tp1 <- tp1*multiplier
  # tp2 <- tp2*multiplier
  # tp3 <- tp3*multiplier
  # tp4 <- tp4*multiplier
  
  # update: add effect on hazard
  # Assume tp4 is unaffected by age, sex and ethnicity
  # lambda_adj <- -log(1 - c(tp1, tp2, tp3)) * multiplier
  # tp_adj <- 1 - exp(-lambda_adj)
  # tp1 <- tp_adj[1]
  # tp2 <- tp_adj[2]
  # tp3 <- tp_adj[3] 
  
  # update: use logit 
  logit_tp <- log(c(tp1, tp2, tp3) / (1 - c(tp1, tp2, tp3)))
  logit_tp_adj <- logit_tp + multiplier
  tp_adj <- exp(logit_tp_adj) / (1 + exp(logit_tp_adj))
  tp1 <- tp_adj[1]
  tp2 <- tp_adj[2]
  tp3 <- tp_adj[3]
  ###
  
  # Generate transition times from exponential distributions
  time_12 <- rexp(n_sim, rate = tp_to_rate(tp12))
  time_23 <- rexp(n_sim, rate = tp_to_rate(tp23))
  time_34 <- rexp(n_sim, rate = tp_to_rate(tp34))
  time_diag_1 <- rexp(n_sim, rate = tp_to_rate(tp1))
  time_diag_2 <- rexp(n_sim, rate = tp_to_rate(tp2))
  time_diag_3 <- rexp(n_sim, rate = tp_to_rate(tp3))
  time_diag_4 <- rexp(n_sim, rate = tp_to_rate(tp4))
  
  # Calculate cumulative progression times
  time_to_stage2 <- time_12
  time_to_stage3 <- time_12 + time_23
  time_to_stage4 <- time_12 + time_23 + time_34
  
  # Initialize counters for each possible starting stage
  stage_counts <- rep(0, 4)
  
  for (i in 1:n_sim) {
    # Determine which stage the patient could have started in
    # based on the diagnosis timing and stage
    
    if (diagnosed_stage == 1) {
      # Diagnosed at stage 1 - must have started at stage 1
      # if (time_diag_1[i] >= month_diagnosis) {
      stage_counts[1] <- stage_counts[1] + 1
      # }
    } else if (diagnosed_stage == 2) {
      # Diagnosed at stage 2 - could have started at stage 1 or 2
      if (time_diag_2[i] >= month_diagnosis) {
        # Started at stage 2
        stage_counts[2] <- stage_counts[2] + 1
      } else {
        # } else if (time_to_stage2[i] <= month_diagnosis && 
        #            time_to_stage2[i] + time_diag_2[i] >= month_diagnosis) {
        # Started at stage 1, progressed to stage 2
        stage_counts[1] <- stage_counts[1] + 1
      }
    } else if (diagnosed_stage == 3) {
      # Diagnosed at stage 3 - could have started at stage 1, 2, or 3
      if (time_diag_3[i] >= month_diagnosis) {
        # Started at stage 3
        stage_counts[3] <- stage_counts[3] + 1
      }
      # } else if (time_to_stage3[i] <= month_diagnosis && 
      #            time_to_stage3[i] + time_diag_3[i] >= month_diagnosis) {
      # Progressed to stage 3 before diagnosis
      else if (time_to_stage2[i] <= month_diagnosis - time_diag_3[i]) {
        # suggest enough time to progress to 2 if start from 1
        # within month_diagnosis - time_diag_3[i]
        # So started at stage 1
        stage_counts[1] <- stage_counts[1] + 1
      } else {
        # Started at stage 2
        stage_counts[2] <- stage_counts[2] + 1
      }
      
    } else if (diagnosed_stage == 4) {
      # Diagnosed at stage 4 - could have started at any stage
      if (time_diag_4[i] >= month_diagnosis) {
        # Started at stage 4
        stage_counts[4] <- stage_counts[4] + 1
      } 
      # else if (time_to_stage4[i] <= month_diagnosis && 
      #            time_to_stage4[i] + time_diag_4[i] >= month_diagnosis) {
      # Progressed to stage 4 before diagnosis
      else if (time_23[i] <= month_diagnosis - time_diag_4[i]) {
        # Suggest enough time to progress to 3 if start at 2
        # within month_diagnosis - time_diag_4[i]
        if (time_to_stage3[i] <= month_diagnosis - time_diag_4[i]) {
          # Started at stage 1
          stage_counts[1] <- stage_counts[1] + 1
        } else {
          # Started at stage 2
          stage_counts[2] <- stage_counts[2] + 1
        }
      } else {
        # Started at stage 3
        stage_counts[3] <- stage_counts[3] + 1
      }
    }
  }

  # Convert counts to probabilities
  total_count <- sum(stage_counts)
  if (total_count == 0) {
    return(rep(0.25, 4))  # Equal probability if no valid transitions
  }
  
  return(stage_counts / total_count)
}


# Function to calculate total log-likelihood for individual patient cohort
calculate_cohort_likelihood <- function(params, patient_data, use_month = TRUE, 
                                        plus_param = FALSE,
                                        n_cycles = 24,
                                        n_sim = 10000,
                                        sigma = 5 # = penalise 466 to total_log_likelihood for 1000 people
                                        ) {
  # patient_data should be a data frame with columns:
  # - month: month of diagnosis (0-24) or 12 for 1 year index
  # - diagnosed_stage: 1, 2, 3, or 4 for stages 1-4
  
  if (plus_param){
    # Input validation
    if (length(params) != 10) return(-Inf)
  } else {
    if (length(params) != 7) return(-Inf)
  }
  
  if (any(is.na(params)) || any(is.null(params))) return(-Inf)
  
  # Extract parameters
  tp12 <- params[1]  # undiagnosed stage 1 -> undiagnosed stage 2
  tp23 <- params[2]  # undiagnosed stage 2 -> undiagnosed stage 3
  tp34 <- params[3]  # undiagnosed stage 3 -> undiagnosed stage 4
  tp1 <- params[4]   # undiagnosed stage 1 -> diagnosed stage 1
  tp2 <- params[5]   # undiagnosed stage 2 -> diagnosed stage 2
  tp3 <- params[6]   # undiagnosed stage 3 -> diagnosed stage 3
  tp4 <- params[7]   # undiagnosed stage 4 -> diagnosed stage 4
  
  # Create transition matrix
  P <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
  
  # Pre-compute simulations for all starting states
  sim_results_origin <- list()
  for (start_state in 1:4) {
    sim_results_origin[[start_state]] <- simulate_markov(P, initial_state = start_state, n_cycles = n_cycles)
  }
  
  if (plus_param){
    
    beta_age70plus <- params[8]
    beta_female <- params[9]
    beta_nonwhite <- params[10]
    
    # give penalty to extreme three beta values
    # the three beta range between negative values to zero
    # penalise large negative values
    n_patient <- nrow(patient_data)
    penalty <- (beta_age70plus + beta_female + beta_nonwhite)^2 * n_patient/(2 * sigma^2)
    
    # Group patients by month and diagnosed stage
    patient_summary <- patient_data %>%
      group_by(month, diagnosed_stage, age70plus, female, nonwhite) %>%
      summarise(count = n(), .groups = 'drop')
  } else {
    # Group patients by month and diagnosed stage
    patient_summary <- patient_data %>%
      group_by(month, diagnosed_stage) %>%
      summarise(count = n(), .groups = 'drop')
    penalty <- 0
  }
  
  # Calculate total log-likelihood
  total_log_likelihood <- 0
  
  for (i in 1:nrow(patient_summary)) {
    month <- patient_summary$month[i]
    diagnosed_stage <- patient_summary$diagnosed_stage[i]
    count <- patient_summary$count[i]
    
    if (month < 0 || month > n_cycles || diagnosed_stage < 1 || diagnosed_stage > 4) next
    
    multiplier <- 0
    
    if (plus_param){
      female <- patient_summary$female[i]
      nonwhite <- patient_summary$nonwhite[i]
      age70plus <- patient_summary$age70plus[i]
      
      # They are multiplicative. 
      # when age70plus/female/nonwhite = 0, the effect = 1
      # multiplier <- (beta_age70plus^age70plus) * (beta_female^female) * (beta_nonwhite^nonwhite)
      
      # update: multiplicative multiplier leads to extreme estimation
      # add effect on hazard
      # multiplier <- exp(beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite)
      
      # update2: multiplicative on hazard still push age and sex to extreme
      # try logit
      multiplier <- beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite
      
      # Create transition matrix
      PP <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4, multiplier = multiplier)
      
      # Pre-compute simulations for all starting states
      sim_results <- list()
      for (start_state in 1:4) {
        sim_results[[start_state]] <- simulate_markov(PP, initial_state = start_state, n_cycles = n_cycles)
      }
    } else {
      sim_results <- sim_results_origin
    }
    
    # Calculate stage distribution probabilities at entry
    stage_probs <- initial_stage_distribution(month, diagnosed_stage, 
                                              tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                              multiplier = multiplier,
                                              n_sim = n_sim)
    
    # Calculate likelihood for this observation
    prob_total <- 0
    
    if (use_month){# Target is the probability of diagnosed at stage X at month Y
      for (start_stage in 1:4) {
        
        sim_data <- sim_results[[start_stage]]
        
        
        if (month + 1 <= nrow(sim_data)) {
          if (month == 0) {
            # At month 0, use cumulative probability of diagnosis at this stage
            prob_scenario <- sim_data[month + 1, diagnosed_stage + 4]  # +4 for diagnosed states
          } else {
            # For month > 0, probability of being undiagnosed at this stage at month-1
            # AND transitioning to diagnosed at this stage at month
            prob_undiag_prev <- sim_data[month, diagnosed_stage]
            trans_prob <- switch(diagnosed_stage, # select the corresponding tp from the four
                                 tp1, tp2, tp3, tp4)
            prob_scenario <- prob_undiag_prev * trans_prob
          }
          prob_total <- prob_total + stage_probs[start_stage] * prob_scenario
        }
      }
    } else {# Target is the probability of diagnosed at stage X overall
      for (start_stage in 1:4) {
        
        sim_data <- sim_results[[start_stage]]
        
        # just look at the row 25, as the numbers are accumulative
        prob_scenario <- sum(sim_data[25, diagnosed_stage + 4])/sum(sim_data[25, 1:8])
        
        prob_total <- prob_total + stage_probs[start_stage] * prob_scenario
        
      }
    }
    
    prob_total <- max(prob_total, 1e-10)
    total_log_likelihood <- total_log_likelihood + count * log(prob_total)
  }
  
  total_log_likelihood <- total_log_likelihood - penalty
  
  return(total_log_likelihood)
}


# Function to save intermediate results
save_intermediate_results <- function(result, filename_base, iteration = NULL) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  if (!is.null(iteration)) {
    filename <- paste0(filename_base, "_iter", iteration, "_", timestamp, ".rds")
  } else {
    filename <- paste0(filename_base, "_", timestamp, ".rds")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists("output")) dir.create("output")
  if (!dir.exists("output/optim")) dir.create("output/optim")
  
  saveRDS(result, file = file.path("output", "optim", filename))
  cat("Results saved to:", filename, "\n")
  return(filename)
}

# Robust parameter validation function
validate_parameters <- function(params) {
  # if (length(params) != 7) return(FALSE)
  if (any(is.na(params)) || any(is.null(params))) return(FALSE)
  
  tp12 <- params[1]
  tp23 <- params[2]
  tp34 <- params[3]
  tp1 <- params[4]
  tp2 <- params[5]
  tp3 <- params[6]
  tp4 <- params[7]
  
  # Check bounds
  if (any(params[1:7] < 0.001) || any(params[1:7] > 1)) return(FALSE)
  
  # Check constraints: sum of outgoing probabilities <= 1
  if ((tp12 + tp1) > 1 || (tp23 + tp2) > 1 || (tp34 + tp3) > 1 || tp4 > 1) {
    return(FALSE)
  }
  
  # Check biological constraint: diagnosis rates should increase with stage
  if (!(tp1 <= tp2 && tp2 <= tp3 && tp3 <= tp4)) {
    return(FALSE)
  }
  
  # Progression rate should accelerate gradually
  if (!(tp12 <= tp23 && tp23 <= tp34)) {
    return(FALSE)
  }
  
  return(TRUE)
}


# Helper function to suggest TPs from dwelling years
dwell2tps <- function(dw_list){
  
  stage1_lower = dw_list["stage1_lower"] 
  stage1_upper = dw_list["stage1_upper"]
  stage2_lower = dw_list["stage2_lower"]
  stage2_upper = dw_list["stage2_upper"]
  stage3_lower = dw_list["stage3_lower"]
  stage3_upper = dw_list["stage3_upper"]
  stage4_lower = dw_list["stage4_lower"]
  stage4_upper =dw_list["stage4_upper"]
  
  # Dwelling times in years for each stage
  # sample dwell years
  # the later stage dwell year is not longer than early stage dwell year. 
  par11 <- runif(1, stage1_lower, stage1_upper)
  par21 <- runif(1, min(stage2_lower, min(par11, stage2_upper)), max(stage2_lower, min(par11, stage2_upper)))
  par31 <- runif(1, min(stage3_lower, min(par21, stage3_upper)), max(stage3_lower, min(par21, stage3_upper)))
  par41 <- runif(1, min(stage4_lower, min(par31, stage4_upper)), max(stage4_lower, min(par31, stage4_upper)))
  
  dwell_years = c(par11, par21, par31, par41)
  
  stage_names <- c("Stage 1", "Stage 2", "Stage 3", "Stage 4")
  
  # Convert to months
  dwell_months <- dwell_years * 12
  
  # Calculate transition rates (hazard rates)
  # For exponential distribution: rate = 1 / mean_time
  lambda <- 1 / dwell_months
  
  # Calculate monthly transition probabilities
  # P(transition in 1 month) = 1 - exp(-lambda * 1)
  p_monthly <- 1 - exp(-lambda * 1)
  
  # Create results dataframe
  rt <- data.frame(
    Stage = stage_names,
    Dwell_Years = dwell_years,
    Dwell_Months = dwell_months,
    Transition_Rate = round(lambda, 4),
    Monthly_Prob = round(p_monthly, 4),
    Monthly_Prob_Pct = round(p_monthly * 100, 2)
  )
 
  return(rt)
}

# p_monthly <- 0.059
# (1/(-log(1 - p_monthly)))/12

# Sample beta for the three additional parameters for age, sex and ethnicity
# 80% from Gaussian, informed by AFT model on diagnostic interval
# 20% from uniform allowing flexibility 
sample_beta <- function(prior_mean, prior_sd, lower, upper, uniform_prob = 0.1) {
  if(runif(1) < uniform_prob) {
    # Uniform exploration
    beta <- runif(1, lower, upper)
  } else {
    # Gaussian sampling centered on AFT estimates
    beta <- rnorm(1, mean = prior_mean, sd = prior_sd)
    # Truncate to bounds
    beta <- pmax(lower, pmin(upper, beta))
  }
  return(beta)
}


# Helper function to suggest good starting parameters
suggest_starting_parameters <- function(stage_prop, 
                                        dw_list,
                                        # args below are for plus age, sex and eth
                                        plus_param = FALSE,
                                        suggest_logit_mean = c(-0.16, -0.19, -0.09),
                                        suggest_logit_sd = c(0.5, 0.5, 0.5),
                                        lower = -1.5,
                                        upper = 0.2
                                        ) {

  suggested_tps <- dwell2tps(dw_list)["Monthly_Prob"]
  
  # Suggest parameters based on data characteristics
  # Progression rates between stages
  suggested_tp12 <- suggested_tps$Monthly_Prob[1]
  suggested_tp23 <- suggested_tps$Monthly_Prob[2]
  suggested_tp34 <- suggested_tps$Monthly_Prob[3]
  
  # Diagnosis rates (should increase with stage)
  s4_diag_rate <- suggested_tps$Monthly_Prob[4]
  suggested_tp4 <- s4_diag_rate
  
  # get distribution of stage at diagnosis
  # use ratio between neighbouring stages * 2 as the upper bound
  # 1 as the lower bound
  down_scale34 <- runif(1, 1, 2*max(1, stage_prop[4]/stage_prop[3]))
  down_scale23 <- runif(1, 1, 2*max(1, stage_prop[3]/stage_prop[2]))
  down_scale12 <- runif(1, 1, 2*max(1, stage_prop[2]/stage_prop[1]))

  suggested_tp3 <- suggested_tp4/down_scale34
  suggested_tp2 <- suggested_tp3/down_scale23
  suggested_tp1 <- suggested_tp2/down_scale12
  
  suggested_params <- c(suggested_tp12, suggested_tp23, suggested_tp34, 
                        suggested_tp1, suggested_tp2, suggested_tp3, suggested_tp4)
  
  if (plus_param){

    # ln_HR_min <- log(HR_min)
    # ln_HR_max <- log(HR_max)
    # 
    # beta_age70plus <- runif(1, ln_HR_min, ln_HR_max)
    # beta_female <- runif(1, ln_HR_min, ln_HR_max)
    # beta_nonwhite <- runif(1, ln_HR_min, ln_HR_max)
    
    beta_age70plus <- sample_beta(
      prior_mean = suggest_logit_mean[1], 
      prior_sd = suggest_logit_sd[1], 
      lower = lower, 
      upper = upper
    )
    
    beta_female <- sample_beta(
      prior_mean = suggest_logit_mean[2], 
      prior_sd = suggest_logit_sd[2], 
      lower = lower, 
      upper = upper
    )
    
    beta_nonwhite <- sample_beta(
      prior_mean = suggest_logit_mean[3], 
      prior_sd = suggest_logit_sd[3], 
      lower = lower, 
      upper = upper
    )
    
    suggested_params_plus <- c(suggested_params, beta_age70plus, beta_female, beta_nonwhite)
    
    return(suggested_params_plus)
  }
  
  return(suggested_params)
}

# Helper function to set reasonable parameter bounds
get_parameter_bounds <- function(plus_param) {
  # Lower bounds (all parameters must be positive)
  lower_bounds <- rep(0.01, 7)
  
  # Upper bounds
  upper_bounds <- c(
    0.5,   # tp12: stage 1->2 progression
    0.5,   # tp23: stage 2->3 progression  
    0.5,   # tp34: stage 3->4 progression
    0.5,   # tp1: stage 1 diagnosis rate
    0.5,   # tp2: stage 2 diagnosis rate
    0.5,   # tp3: stage 3 diagnosis rate
    0.5    # tp4: stage 4 diagnosis rate
  )
  
  if (plus_param){
    lower_bounds <- c(lower_bounds, c(-1.5, -1.5, -1.5))
    upper_bounds <- c(upper_bounds, c(0.2, 0.2, 0.2))
  }
  
  return(list(lower = lower_bounds, upper = upper_bounds))
}

# Function to convert transition probability matrix to rate matrix
# The simple approximation Q ≈ (P - I)/Δt works well when transition probabilities are small (< 0.2)
# in effect, when tp is small, tp ≈ rate
probability_to_rate_matrix <- function(P, time_step = 1, simple = TRUE) {
  # P: transition probability matrix
  # time_step: time interval (default = 1)
  # This works well when transition probabilities are small
  I <- diag(nrow(P))
  Q <- (P - I) / time_step
  
  # Function using matrix logarithm (more accurate for larger probabilities)
  if(!simple){
    # Requires expm package for matrix operations
    if (!requireNamespace("expm", quietly = TRUE)) {
      stop("Package 'expm' is required for exact matrix logarithm calculation")
    }
    
    # Q = logm(P) / time_step
    Q <- expm::logm(P) / time_step
  }
  
  return(Q)
}

# Helper function to generate valid random parameters
generate_valid_params <- function(lower_bounds, upper_bounds, max_attempts = 10000) {
  
  # Create a function for 3-decimal random sampling
  runif3 <- function(n, min, max) {
    grid <- seq(min, max, by = 0.005)
    sample(grid, n, replace = TRUE)
  }
  
  for (attempt in 1:max_attempts) {
    # Generate parameters that respect ordering constraints
    
    # Progression rates (tp12 <= tp23 <= tp34)
    tp12 <- runif3(1, lower_bounds[1], upper_bounds[1])
    tp23 <- runif3(1, tp12, upper_bounds[2])  # Must be >= tp12
    tp34 <- runif3(1, tp23, upper_bounds[3])  # Must be >= tp23
    
    # Diagnosis rates (tp1 <= tp2 <= tp3 <= tp4)
    tp1 <- runif3(1, lower_bounds[4], upper_bounds[4])
    tp2 <- runif3(1, tp1, upper_bounds[5])    # Must be >= tp1
    tp3 <- runif3(1, tp2, upper_bounds[6])    # Must be >= tp2
    tp4 <- runif3(1, tp3, upper_bounds[7])    # Must be >= tp3
    
    params <- c(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
    
    # Check sum constraints: outgoing probabilities <= 1
    if ((tp12 + tp1) <= 1 && (tp23 + tp2) <= 1 && 
        (tp34 + tp3) <= 1 && tp4 <= 1) {
      return(params)
    }
  }
  
  # If max_attempts reached, return NULL
  warning("Could not generate valid parameters after ", max_attempts, " attempts")
  return(NULL)
}



tp_to_rate <- function(tp){
  
  lambda <- -log(1 - tp)
  
  return(lambda)
}

