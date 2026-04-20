# Parameters for cancer natural history model
# Cancer type: galb, oeso, panc, stom

# Optimisation ------------------------------------------------------------

# Functions of optimisation of natural history parameters - 4-stage version

library(tidyverse)
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

create_transition_matrix <- function(tp12, tp23, tp34, tp1, tp2, tp3, tp4) {
  # tp12: undiagnosed stage 1 -> undiagnosed stage 2
  # tp23: undiagnosed stage 2 -> undiagnosed stage 3
  # tp34: undiagnosed stage 3 -> undiagnosed stage 4
  # tp1: undiagnosed stage 1 -> diagnosed stage 1
  # tp2: undiagnosed stage 2 -> diagnosed stage 2
  # tp3: undiagnosed stage 3 -> diagnosed stage 3
  # tp4: undiagnosed stage 4 -> diagnosed stage 4
  
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
calculate_stage_distribution <- function(month_diagnosis, diagnosed_stage, 
                                         tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                         n_sim = 10000, seed_num = 1234) {
  set.seed(seed_num)
  
  # Generate transition times from exponential distributions
  time_12 <- rexp(n_sim, rate = tp12)
  time_23 <- rexp(n_sim, rate = tp23)
  time_34 <- rexp(n_sim, rate = tp34)
  time_diag_1 <- rexp(n_sim, rate = tp1)
  time_diag_2 <- rexp(n_sim, rate = tp2)
  time_diag_3 <- rexp(n_sim, rate = tp3)
  time_diag_4 <- rexp(n_sim, rate = tp4)
  
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
        # Started at stage 1
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
        if (time_to_stage2[i] <= month_diagnosis - time_diag_4[i]) {
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
calculate_cohort_likelihood <- function(params, patient_data) {
  # patient_data should be a data frame with columns:
  # - month: month of diagnosis (0-24)
  # - diagnosed_stage: 1, 2, 3, or 4 for stages 1-4
  
  # Input validation
  if (length(params) != 7) return(-Inf)
  
  if (any(is.na(params)) || any(is.null(params))) return(-Inf)
  
  # Extract parameters
  tp12 <- params[1]  # undiagnosed stage 1 -> undiagnosed stage 2
  tp23 <- params[2]  # undiagnosed stage 2 -> undiagnosed stage 3
  tp34 <- params[3]  # undiagnosed stage 3 -> undiagnosed stage 4
  tp1 <- params[4]   # undiagnosed stage 1 -> diagnosed stage 1
  tp2 <- params[5]   # undiagnosed stage 2 -> diagnosed stage 2
  tp3 <- params[6]   # undiagnosed stage 3 -> diagnosed stage 3
  tp4 <- params[7]   # undiagnosed stage 4 -> diagnosed stage 4
  
  # Parameter constraints
  if (any(params < 0) || any(params > 1)) return(-Inf)
  
  # Constraint: sum of outgoing probabilities <= 1 for each state
  if ((tp12 + tp1) > 1 || (tp23 + tp2) > 1 || (tp34 + tp3) > 1 || tp4 > 1) {
    return(-Inf)
  }
  
  # Biological constraint: diagnosis rates should generally increase with stage
  if (!(tp1 <= tp2 && tp2 <= tp3 && tp3 <= tp4)) {
    return(-Inf)
  }
  
  # Create transition matrix
  P <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
  
  # Pre-compute simulations for all starting states
  sim_results <- list()
  for (start_state in 1:4) {
    sim_results[[start_state]] <- simulate_markov(P, initial_state = start_state, n_cycles = 24)
  }
  
  # Group patients by month and diagnosed stage
  patient_summary <- patient_data %>%
    group_by(month, diagnosed_stage) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate total log-likelihood
  total_log_likelihood <- 0
  
  for (i in 1:nrow(patient_summary)) {
    month <- patient_summary$month[i]
    diagnosed_stage <- patient_summary$diagnosed_stage[i]
    count <- patient_summary$count[i]
    
    if (month < 0 || month > 24 || diagnosed_stage < 1 || diagnosed_stage > 4) next
    
    # Calculate stage distribution probabilities at entry
    stage_probs <- calculate_stage_distribution(month, diagnosed_stage, 
                                                tp12, tp23, tp34, tp1, tp2, tp3, tp4)
    
    # Calculate likelihood for this observation
    prob_total <- 0
    
    for (start_stage in 1:4) {
      if (stage_probs[start_stage] > 0) {
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
    }
    
    prob_total <- max(prob_total, 1e-10)
    total_log_likelihood <- total_log_likelihood + count * log(prob_total)
  }
  
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
  if (length(params) != 7) return(FALSE)
  if (any(is.na(params)) || any(is.null(params))) return(FALSE)
  
  tp12 <- params[1]
  tp23 <- params[2]
  tp34 <- params[3]
  tp1 <- params[4]
  tp2 <- params[5]
  tp3 <- params[6]
  tp4 <- params[7]
  
  # Check bounds
  if (any(params < 0) || any(params > 1)) return(FALSE)
  
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


# Maximum likelihood estimation for individual patient data
estimate_parameters <- function(patient_data, initial_params, 
                                lower_params, upper_params,
                                save_interval = NA,
                                filename_base = "mle_results_4stage",
                                save_final = TRUE) {
  iteration_count <- 0
  
  # Validate initial parameters
  if (!validate_parameters(initial_params)) {
    stop("Initial parameters are invalid. Please check constraints.")
  }
  
  # Test likelihood function with initial parameters
  test_likelihood <- calculate_cohort_likelihood(initial_params, patient_data)
  if (is.infinite(test_likelihood)) {
    warning("Initial parameters produce infinite likelihood. Consider adjusting starting values.")
  }
  
  # Objective function (negative log-likelihood) with progress tracking
  objective <- function(params) {
    iteration_count <<- iteration_count + 1
    
    # Validate parameters first
    if (!validate_parameters(params)) {
      return(1e6)  # Return large but finite penalty
    }
    
    # Try to calculate likelihood with error handling
    likelihood <- tryCatch({
      calculate_cohort_likelihood(params, patient_data)
    }, error = function(e) {
      cat("Error in likelihood calculation at iteration", iteration_count, ":\n")
      cat("Parameters:", params, "\n")
      cat("Error:", e$message, "\n")
      return(-Inf)
    })
    
    # Handle infinite likelihood
    if (is.infinite(likelihood)) {
      return(1e6)  # Return large but finite penalty
    }
    
    # Handle NaN or NA
    if (is.na(likelihood)) {
      return(1e6)  # Return large but finite penalty
    }
    
    # Save intermediate results periodically
    if (!is.na(save_interval) && iteration_count %% save_interval == 0) {
      current_result <- list(
        iteration = iteration_count,
        params = params,
        likelihood = likelihood,
        timestamp = Sys.time()
      )
      
      save_intermediate_results(current_result, 
                                paste0(filename_base, "_intermediate"), 
                                iteration_count)
      
      cat("Iteration", iteration_count, "- Likelihood:", likelihood, 
          "- Params:", round(params, 4), "\n")
    }
    
    return(-likelihood)
  }
  
  cat("Starting optimization for 4-stage model...\n")
  cat("Initial parameters:", initial_params, "\n")
  cat("Parameter bounds: [", lower_params, "] to [", upper_params, "]\n")
  
  start_time <- Sys.time()
  
  # Optimization with constraints
  result <- tryCatch({
    optim(
      par = initial_params,
      fn = objective,
      method = "L-BFGS-B",
      lower = lower_params,
      upper = upper_params,
      control = list(
        trace = 1, 
        maxit = 3000,  # Increased for more complex model
        factr = 1e12,
        pgtol = 1e-8,
        ndeps = rep(1e-8, length(initial_params))
      )
    )
  }, error = function(e) {
    cat("Optimization failed with error:", e$message, "\n")
    
    # Try alternative optimization method
    cat("Trying alternative optimization method (Nelder-Mead)...\n")
    optim(
      par = initial_params,
      fn = objective,
      method = "Nelder-Mead",
      control = list(trace = 1, maxit = 2000)
    )
  })
  
  end_time <- Sys.time()
  
  # Validate final results
  if (!validate_parameters(result$par)) {
    warning("Final parameters violate constraints. Results may not be reliable.")
  }
  
  # Save final results
  final_result <- list(
    mle_result = result,
    total_iterations = iteration_count,
    computation_time = end_time - start_time,
    final_params = result$par,
    final_likelihood = -result$value,
    convergence_code = result$convergence,
    convergence_message = switch(as.character(result$convergence),
                                 "0" = "Successful convergence",
                                 "1" = "Maximum iterations reached",
                                 "51" = "Warning from L-BFGS-B",
                                 "52" = "Error from L-BFGS-B",
                                 paste("Unknown convergence code:", result$convergence)),
    parameter_names = c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
  )
  
  final_filename <- NULL
  if (save_final) {
    final_filename <- tryCatch({
      save_intermediate_results(final_result, filename_base)
    }, error = function(e) {
      cat("Warning: Could not save final results:", e$message, "\n")
      return(NULL)
    })
  }
  
  cat("\nOptimization complete!\n")
  cat("Convergence:", final_result$convergence_message, "\n")
  cat("Total iterations:", iteration_count, "\n")
  cat("Total time:", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")
  
  cat("Final parameters:\n")
  for (i in 1:length(result$par)) {
    cat(sprintf("  %s: %.4f\n", final_result$parameter_names[i], result$par[i]))
  }
  cat("Final likelihood:", round(-result$value, 4), "\n")
  
  if (!is.null(final_filename)) {
    cat("Final results saved to:", final_filename, "\n")
  }
  
  return(final_result)
}

# Helper function to suggest good starting parameters
suggest_starting_parameters <- function(patient_data) {
  # Basic analysis of the data
  stage_counts <- table(patient_data$diagnosed_stage)
  total_count <- nrow(patient_data)
  
  # Calculate proportion diagnosed at each stage
  stage_props <- stage_counts / total_count
  
  cat("Data summary:\n")
  cat("Total patients:", total_count, "\n")
  for (i in 1:4) {
    if (i %in% names(stage_counts)) {
      cat(sprintf("Stage %d: %d (%.1f%%)\n", i, stage_counts[as.character(i)], 
                  stage_props[as.character(i)] * 100))
    } else {
      cat(sprintf("Stage %d: 0 (0.0%%)\n", i))
    }
  }
  
  # Suggest parameters based on data characteristics
  # Progression rates between stages
  suggested_tp12 <- 0.025
  suggested_tp23 <- 0.05
  suggested_tp34 <- 0.075
  
  # Diagnosis rates (should increase with stage)
  base_diag_rate <- 0.02
  suggested_tp1 <- base_diag_rate
  suggested_tp2 <- base_diag_rate * 1.5
  suggested_tp3 <- base_diag_rate * 2.5
  suggested_tp4 <- base_diag_rate * 4.0
  
  # Adjust based on observed stage distribution
  if (length(stage_props) > 0) {
    # If mostly late-stage diagnoses, increase later diagnosis rates
    if (sum(stage_props[c("3", "4")], na.rm = TRUE) > 0.6) {
      suggested_tp3 <- suggested_tp3 * 1.5
      suggested_tp4 <- suggested_tp4 * 1.5
    }
  }
  
  suggested_params <- c(suggested_tp12, suggested_tp23, suggested_tp34, 
                        suggested_tp1, suggested_tp2, suggested_tp3, suggested_tp4)
  
  cat("\nSuggested starting parameters:\n")
  param_names <- c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
  for (i in 1:length(suggested_params)) {
    cat(sprintf("  %s: %.4f\n", param_names[i], suggested_params[i]))
  }
  
  return(suggested_params)
}

# Helper function to set reasonable parameter bounds
get_parameter_bounds <- function() {
  # Lower bounds (all parameters must be positive)
  lower_bounds <- rep(0.001, 7)
  
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
  
  return(list(lower = lower_bounds, upper = upper_bounds))
}



# Validation --------------------------------------------------------------

# Function to estimate starting stage distribution for each cancer site
estimate_starting_distribution <- function(d4s, params) {
  # Group by cancer site and calculate stage distribution for each
  sites <- unique(d4s$site)
  starting_distributions <- list()
  
  for (st in sites) {
    site_data <- d4s[d4s$site == st, ]
    
    optimized_params <- as.numeric(params[params$site==st, grepl("tp", colnames(params))])
    
    # Calculate weighted average of starting distributions across all patients
    stage_dist_matrix <- matrix(0, nrow = nrow(site_data), ncol = 4)
    
    for (i in 1:nrow(site_data)) {
      stage_dist_matrix[i, ] <- calculate_stage_distribution(
        month_diagnosis = site_data$month[i],
        diagnosed_stage = site_data$diagnosed_stage[i],
        tp12 = optimized_params[1], tp23 = optimized_params[2], tp34 = optimized_params[3],
        tp1 = optimized_params[4], tp2 = optimized_params[5], 
        tp3 = optimized_params[6], tp4 = optimized_params[7]
      )
    }
    
    # Average across patients for this cancer site
    starting_distributions[[st]] <- colMeans(stage_dist_matrix)
  }
  
  return(starting_distributions)
}

# Microsimulation function for single patient
simulate_single_patient <- function(P, starting_stage_probs) {
  # Sample starting stage based on probabilities
  start_stage <- sample(1:4, 1, prob = starting_stage_probs)
  
  current_state <- start_stage
  
  for (cycle in 1:24) {
    # If already diagnosed, return diagnosis stage
    if (current_state > 4) {
      return(current_state - 4)  # Convert to stage (5->1, 6->2, 7->3, 8->4)
    }
    
    # Sample next state based on transition probabilities
    next_state <- sample(1:8, 1, prob = P[current_state, ])
    current_state <- next_state
    
    # If transitioned to diagnosed state, return diagnosis stage
    if (current_state > 4) {
      return(current_state - 4)
    }
  }
  
  # If still undiagnosed after 24 cycles, return current stage
  return(current_state)
}

# Main prediction function
predict_stage_distribution <- function(patient_data, params, n_sim = 1000) {
  
  # Estimate starting distributions by cancer site
  starting_distributions <- estimate_starting_distribution(patient_data, params)
  
  # Initialize results storage
  sites <- names(starting_distributions)
  predicted_distributions <- list()
  
  for (st in sites) {
    site_results <- matrix(0, nrow = n_sim, ncol = 4)
    site_data <- patient_data[patient_data$site == st, ]
    n_patients <- nrow(site_data)
    
    # Create transition matrix
    optimized_params <- as.numeric(params[params$site==st, grepl("tp", colnames(params))])
    P <- create_transition_matrix(
      optimized_params[1], optimized_params[2], optimized_params[3],
      optimized_params[4], optimized_params[5], optimized_params[6], optimized_params[7]
    )
    
    cat("Simulating", st, "cancer (", n_patients, "patients)...\n")
    
    for (sim in 1:n_sim) {
      stage_counts <- rep(0, 4)
      
      # Simulate each patient in this cancer site
      for (patient in 1:n_patients) {
        diagnosed_stage <- simulate_single_patient(P, starting_distributions[[st]])
        stage_counts[diagnosed_stage] <- stage_counts[diagnosed_stage] + 1
      }
      
      # Store proportions for this simulation
      site_results[sim, ] <- stage_counts / n_patients
    }
    
    # Calculate mean and confidence intervals across simulations
    predicted_distributions[[st]] <- list(
      mean = colMeans(site_results),
      ci_lower = apply(site_results, 2, quantile, 0.025),
      ci_upper = apply(site_results, 2, quantile, 0.975),
      observed = as.numeric(table(factor(site_data$diagnosed_stage, levels = 1:4)) / nrow(site_data))
    )
  }
  
  return(predicted_distributions)
}


# Function to compare predicted vs observed distributions
compare_distributions <- function(prediction_results) {
  comparison_df <- data.frame()
  
  for (site in names(prediction_results)) {
    site_data <- prediction_results[[site]]
    
    for (stage in 1:4) {
      comparison_df <- rbind(comparison_df, data.frame(
        cancer_site = site,
        stage = stage,
        predicted = site_data$mean[stage],
        predicted_lower = site_data$ci_lower[stage],
        predicted_upper = site_data$ci_upper[stage],
        observed = site_data$observed[stage],
        difference = site_data$mean[stage] - site_data$observed[stage]
      ))
    }
  }
  
  return(comparison_df)
}

# Visualization function
plot_comparison <- function(comparison_df) {
  library(ggplot2)
  
  p <- ggplot(comparison_df, aes(x = factor(stage))) +
    geom_col(aes(y = observed, fill = "Observed"), alpha = 0.7, position = "dodge") +
    geom_col(aes(y = predicted, fill = "Predicted"), alpha = 0.7, position = "dodge2") +
    geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper), 
                  width = 0.2, position = position_dodge2(width = 0.9)) +
    facet_wrap(~cancer_site, scales = "free") +
    labs(x = "Cancer Stage", y = "Proportion", 
         title = "Predicted vs Observed Stage Distribution by Cancer Site") +
    scale_fill_manual(values = c("Observed" = "steelblue", "Predicted" = "coral")) +
    theme_minimal()
  
  return(p)
}


# Calibration -------------------------------------------------------------

# Calibration Functions for Cancer Natural History Model Parameters

library(tidyverse)

# Enhanced calibration adjustment factors calculation
calculate_calibration_factors <- function(prediction_results, target_weight = 0.3) {
  
  # Calculate comprehensive adjustment factors for both progression and diagnosis rates
  # to better match observed stage distributions
  # 
  # Args:
  #   prediction_results: Output from predict_stage_distribution() function
  #   target_weight: Weight given to adjustment (0-1, lower = more conservative)
  # 
  # Returns:
  #   List of adjustment factors for each cancer site with both diagnosis and progression adjustments
  
  calibration_factors <- list()
  
  for (site in names(prediction_results)) {
    site_data <- prediction_results[[site]]
    
    # Calculate ratio of observed to predicted for each stage
    obs_pred_ratios <- site_data$observed / site_data$mean
    
    # Handle division by zero or very small predicted values
    obs_pred_ratios[site_data$mean < 0.001] <- 1.0
    obs_pred_ratios[is.infinite(obs_pred_ratios)] <- 1.0
    obs_pred_ratios[is.na(obs_pred_ratios)] <- 1.0
    
    # Calculate stage-specific adjustments
    # If early stages are over-predicted, reduce their diagnosis rates and increase progression
    # If late stages are under-predicted, increase their diagnosis rates
    
    diagnosis_adjustments <- rep(1.0, 4)
    progression_adjustments <- rep(1.0, 3)  # tp12, tp23, tp34
    
    # Diagnosis rate adjustments (direct from obs/pred ratios)
    diagnosis_adjustments <- 1 + target_weight * (obs_pred_ratios - 1)
    
    # Progression rate adjustments (inverse relationship with early stage over-prediction)
    if (obs_pred_ratios[1] < 1) {  # Stage 1 over-predicted
      progression_adjustments[1] <- 1 + target_weight * (1/obs_pred_ratios[1] - 1)  # Increase tp12
    }
    if (obs_pred_ratios[2] < 1) {  # Stage 2 over-predicted
      progression_adjustments[2] <- 1 + target_weight * (1/obs_pred_ratios[2] - 1)  # Increase tp23
    }
    if (obs_pred_ratios[3] < 1) {  # Stage 3 over-predicted
      progression_adjustments[3] <- 1 + target_weight * (1/obs_pred_ratios[3] - 1)  # Increase tp34
    }
    
    # Additional logic: if later stages are under-predicted, reduce progression to earlier stages
    if (obs_pred_ratios[4] > 1) {  # Stage 4 under-predicted
      progression_adjustments[1] <- progression_adjustments[1] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.3)
      progression_adjustments[2] <- progression_adjustments[2] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.5)
      progression_adjustments[3] <- progression_adjustments[3] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.7)
    }
    
    # Ensure factors are reasonable (between 0.3 and 3.0)
    diagnosis_adjustments <- pmax(0.3, pmin(3.0, diagnosis_adjustments))
    progression_adjustments <- pmax(0.3, pmin(3.0, progression_adjustments))
    
    calibration_factors[[site]] <- list(
      diagnosis = diagnosis_adjustments,
      progression = progression_adjustments
    )
    
    cat("Calibration factors for", site, ":\n")
    cat("  Diagnosis rates (tp1-tp4):\n")
    for (stage in 1:4) {
      cat(sprintf("    Stage %d: %.3f (obs: %.3f, pred: %.3f)\n", 
                  stage, diagnosis_adjustments[stage], 
                  site_data$observed[stage], site_data$mean[stage]))
    }
    cat("  Progression rates (tp12, tp23, tp34):\n")
    prog_names <- c("tp12", "tp23", "tp34")
    for (i in 1:3) {
      cat(sprintf("    %s: %.3f\n", prog_names[i], progression_adjustments[i]))
    }
    cat("\n")
  }
  
  return(calibration_factors)
}


# Enhanced parameter calibration with comprehensive constraint checking
calibrate_parameters <- function(original_params, calibration_factors) {
  
  # Apply calibration adjustments to both progression and diagnosis rate parameters
  # with comprehensive constraint validation
  # 
  # Args:
  #   original_params: DataFrame with optimized parameters for each site
  #   calibration_factors: Output from calculate_calibration_factors()
  # 
  # Returns:
  #   DataFrame with calibrated parameters
  
  calibrated_params <- original_params
  
  for (site in names(calibration_factors)) {
    if (site %in% calibrated_params$site) {
      site_row <- which(calibrated_params$site == site)
      diagnosis_factors <- calibration_factors[[site]]$diagnosis
      progression_factors <- calibration_factors[[site]]$progression
      
      # Store original values for comparison
      original_values <- as.numeric(calibrated_params[site_row, c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")])
      
      # Apply factors to progression rate parameters (tp12, tp23, tp34)
      calibrated_params[site_row, "tp12"] <- calibrated_params[site_row, "tp12"] * progression_factors[1]
      calibrated_params[site_row, "tp23"] <- calibrated_params[site_row, "tp23"] * progression_factors[2]
      calibrated_params[site_row, "tp34"] <- calibrated_params[site_row, "tp34"] * progression_factors[3]
      
      # Apply factors to diagnosis rate parameters (tp1, tp2, tp3, tp4)
      calibrated_params[site_row, "tp1"] <- calibrated_params[site_row, "tp1"] * diagnosis_factors[1]
      calibrated_params[site_row, "tp2"] <- calibrated_params[site_row, "tp2"] * diagnosis_factors[2]
      calibrated_params[site_row, "tp3"] <- calibrated_params[site_row, "tp3"] * diagnosis_factors[3]
      calibrated_params[site_row, "tp4"] <- calibrated_params[site_row, "tp4"] * diagnosis_factors[4]
      
      # Ensure all parameters are within (0, 1) bounds
      param_cols <- c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
      for (col in param_cols) {
        calibrated_params[site_row, col] <- pmax(0.01, pmin(0.999, calibrated_params[site_row, col]))
      }
      
      # Apply constraints from validate_parameters function
      calibrated_params <- apply_parameter_constraints(calibrated_params, site_row)
      
      # Report changes
      cat("Calibrated parameters for", site, ":\n")
      new_values <- as.numeric(calibrated_params[site_row, param_cols])
      param_names <- c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
      
      for (i in 1:length(param_names)) {
        old_val <- original_values[i]
        new_val <- new_values[i]
        change_pct <- if(old_val != 0) 100 * (new_val - old_val) / old_val else 0
        cat(sprintf("  %s: %.4f -> %.4f (change: %+.1f%%)\n", 
                    param_names[i], old_val, new_val, change_pct))
      }
      cat("\n")
    }
  }
  
  return(calibrated_params)
}


# Function to apply parameter constraints as specified in validate_parameters
apply_parameter_constraints <- function(params_df, site_row) {
  
  # Extract current parameter values
  tp12 <- params_df[site_row, "tp12"]
  tp23 <- params_df[site_row, "tp23"]
  tp34 <- params_df[site_row, "tp34"]
  tp1 <- params_df[site_row, "tp1"]
  tp2 <- params_df[site_row, "tp2"]
  tp3 <- params_df[site_row, "tp3"]
  tp4 <- params_df[site_row, "tp4"]
  
  # Constraint 1: Sum of outgoing probabilities <= 1 for each state
  constraints_violated <- TRUE
  max_iterations <- 50
  iteration <- 0
  
  while (constraints_violated && iteration < max_iterations) {
    iteration <- iteration + 1
    constraints_violated <- FALSE
    
    # Stage 1: tp12 + tp1 <= 1
    if ((tp12 + tp1) > 1) {
      scale_factor <- 0.99 / (tp12 + tp1)
      tp12 <- tp12 * scale_factor
      tp1 <- tp1 * scale_factor
      constraints_violated <- TRUE
    }
    
    # Stage 2: tp23 + tp2 <= 1
    if ((tp23 + tp2) > 1) {
      scale_factor <- 0.99 / (tp23 + tp2)
      tp23 <- tp23 * scale_factor
      tp2 <- tp2 * scale_factor
      constraints_violated <- TRUE
    }
    
    # Stage 3: tp34 + tp3 <= 1
    if ((tp34 + tp3) > 1) {
      scale_factor <- 0.99 / (tp34 + tp3)
      tp34 <- tp34 * scale_factor
      tp3 <- tp3 * scale_factor
      constraints_violated <- TRUE
    }
    
    # Stage 4: tp4 <= 1 (already handled by bounds)
    if (tp4 > 1) {
      tp4 <- 0.99
      constraints_violated <- TRUE
    }
  }
  
  # Constraint 2: Biological constraint - diagnosis rates should increase with stage
  # Adjust gradually to maintain ordering while respecting calibration intent
  diagnosis_rates <- as.numeric(c(tp1, tp2, tp3, tp4))
  
  # If ordering is violated, apply gentle correction
  for (stage in 1:3) {
    if (diagnosis_rates[stage] > diagnosis_rates[stage + 1]) {
      # Average the two rates and add small increment for higher stage
      avg_rate <- (diagnosis_rates[stage] + diagnosis_rates[stage + 1]) / 2
      diagnosis_rates[stage] <- avg_rate * 0.95
      diagnosis_rates[stage + 1] <- avg_rate * 1.05
    }
  }
  
  tp1 <- diagnosis_rates[1]
  tp2 <- diagnosis_rates[2]
  tp3 <- diagnosis_rates[3]
  tp4 <- diagnosis_rates[4]
  
  # Constraint 3: Progression rates should generally increase (optional, gentle enforcement)
  progression_rates <- as.numeric(c(tp12, tp23, tp34))
  for (stage in 1:2) {
    if (progression_rates[stage] > progression_rates[stage + 1] * 1.5) {
      # Only adjust if the violation is significant
      avg_rate <- (progression_rates[stage] + progression_rates[stage + 1]) / 2
      progression_rates[stage] <- avg_rate * 1.1
      progression_rates[stage + 1] <- avg_rate * 1.2
    }
  }
  
  tp12 <- progression_rates[1]
  tp23 <- progression_rates[2]
  tp34 <- progression_rates[3]
  
  # Final bounds check
  tp12 <- pmax(0.001, pmin(0.999, tp12))
  tp23 <- pmax(0.001, pmin(0.999, tp23))
  tp34 <- pmax(0.001, pmin(0.999, tp34))
  tp1 <- pmax(0.001, pmin(0.999, tp1))
  tp2 <- pmax(0.001, pmin(0.999, tp2))
  tp3 <- pmax(0.001, pmin(0.999, tp3))
  tp4 <- pmax(0.001, pmin(0.999, tp4))
  
  # Update the dataframe
  params_df[site_row, "tp12"] <- tp12
  params_df[site_row, "tp23"] <- tp23
  params_df[site_row, "tp34"] <- tp34
  params_df[site_row, "tp1"] <- tp1
  params_df[site_row, "tp2"] <- tp2
  params_df[site_row, "tp3"] <- tp3
  params_df[site_row, "tp4"] <- tp4
  
  # Validate final parameters
  if (!validate_single_site_parameters(tp12, tp23, tp34, tp1, tp2, tp3, tp4)) {
    warning(paste("Parameters for site in row", site_row, "still violate constraints after adjustment"))
  }
  
  return(params_df)
}


# Helper function to validate parameters for a single site
validate_single_site_parameters <- function(tp12, tp23, tp34, tp1, tp2, tp3, tp4) {
  
  # Check bounds
  params <- c(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
  if (any(params < 0) || any(params > 1)) return(FALSE)
  
  # Check constraints: sum of outgoing probabilities <= 1
  if ((tp12 + tp1) > 1.001 || (tp23 + tp2) > 1.001 || (tp34 + tp3) > 1.001 || tp4 > 1.001) {
    return(FALSE)
  }
  
  # Check biological constraint: diagnosis rates should increase with stage
  if (!(tp1 <= tp2 * 1.1 && tp2 <= tp3 * 1.1 && tp3 <= tp4 * 1.1)) {
    return(FALSE)
  }
  
  return(TRUE)
}


# Enhanced iterative calibration with improved convergence criteria
iterative_calibration <- function(patient_data, initial_params, 
                                  max_iterations = 8, 
                                  convergence_threshold = 0.015,
                                  adjustment_weight = 0.25,
                                  n_sim = 1000,
                                  adaptive_weight = TRUE) {
  
  # Perform iterative calibration with enhanced features
  # 
  # Args:
  #   patient_data: Original patient data
  #   initial_params: Starting parameters (output from optimization)
  #   max_iterations: Maximum number of calibration iterations
  #   convergence_threshold: Stop when max absolute difference < this value
  #   adjustment_weight: Initial weight for calibration adjustments (0-1)
  #   n_sim: Number of simulations for prediction
  #   adaptive_weight: Whether to adaptively adjust weight based on progress
  # 
  # Returns:
  #   List containing final calibrated parameters and convergence info
  
  current_params <- initial_params
  iteration <- 0
  convergence_history <- list()
  current_weight <- adjustment_weight
  
  cat("Starting enhanced iterative calibration...\n")
  cat("Initial adjustment weight:", current_weight, "\n\n")
  
  # Get baseline metrics
  baseline_predictions <- predict_stage_distribution(patient_data, current_params, n_sim = n_sim)
  baseline_rmse <- calculate_overall_rmse(baseline_predictions)
  
  repeat {
    iteration <- iteration + 1
    cat("=== Calibration Iteration", iteration, "===\n")
    
    # Predict with current parameters
    predictions <- predict_stage_distribution(patient_data, current_params, n_sim = n_sim)
    
    # Calculate comprehensive metrics
    metrics <- calculate_calibration_metrics(predictions)
    
    cat("Maximum absolute difference:", round(metrics$max_diff, 4), "\n")
    cat("Root Mean Square Error:", round(metrics$rmse, 4), "\n")
    cat("Mean absolute error:", round(metrics$mae, 4), "\n")
    
    # Store convergence metrics
    convergence_history[[iteration]] <- list(
      iteration = iteration,
      max_diff = metrics$max_diff,
      rmse = metrics$rmse,
      mae = metrics$mae,
      adjustment_weight = current_weight,
      params = current_params
    )
    
    # Check convergence
    if (metrics$max_diff < convergence_threshold) {
      cat("Convergence achieved! Max difference < threshold (", convergence_threshold, ")\n")
      break
    }
    
    if (iteration >= max_iterations) {
      cat("Maximum iterations reached (", max_iterations, ")\n")
      break
    }
    
    # Adaptive weight adjustment
    if (adaptive_weight && iteration > 1) {
      prev_rmse <- convergence_history[[iteration-1]]$rmse
      if (metrics$rmse > prev_rmse * 0.95) {  # If not improving significantly
        current_weight <- current_weight * 0.8  # Reduce adjustment weight
        cat("Reducing adjustment weight to:", round(current_weight, 3), "\n")
      }
    }
    
    # Calculate and apply calibration
    cal_factors <- calculate_calibration_factors(predictions, target_weight = current_weight)
    current_params <- calibrate_parameters(current_params, cal_factors)
    
    cat("Moving to next iteration...\n\n")
  }
  
  # Final validation
  cat("=== Final Calibration Results ===\n")
  final_predictions <- predict_stage_distribution(patient_data, current_params, n_sim = n_sim)
  final_comparison <- compare_distributions(final_predictions)
  final_metrics <- calculate_calibration_metrics(final_predictions)
  
  cat("Final metrics:\n")
  cat("  RMSE improvement:", round(baseline_rmse - final_metrics$rmse, 4), "\n")
  cat("  Final max absolute difference:", round(final_metrics$max_diff, 4), "\n")
  cat("  Final RMSE:", round(final_metrics$rmse, 4), "\n")
  
  cat("\nFinal comparison summary:\n")
  print(final_comparison %>% 
          group_by(cancer_site) %>% 
          summarise(
            mean_abs_diff = mean(abs(difference)),
            max_abs_diff = max(abs(difference)),
            rmse = sqrt(mean(difference^2)),
            .groups = 'drop'
          ))
  
  return(list(
    calibrated_params = current_params,
    final_predictions = final_predictions,
    final_comparison = final_comparison,
    convergence_history = convergence_history,
    n_iterations = iteration,
    rmse_improvement = baseline_rmse - final_metrics$rmse
  ))
}


# Helper function to calculate comprehensive calibration metrics
calculate_calibration_metrics <- function(predictions) {
  
  max_diff <- 0
  total_sse <- 0
  total_abs_error <- 0
  n_comparisons <- 0
  
  for (site in names(predictions)) {
    site_data <- predictions[[site]]
    differences <- site_data$mean - site_data$observed
    abs_differences <- abs(differences)
    
    max_diff <- max(max_diff, max(abs_differences))
    total_sse <- total_sse + sum(differences^2)
    total_abs_error <- total_abs_error + sum(abs_differences)
    n_comparisons <- n_comparisons + length(differences)
  }
  
  return(list(
    max_diff = max_diff,
    rmse = sqrt(total_sse / n_comparisons),
    mae = total_abs_error / n_comparisons
  ))
}


# Helper function to calculate overall RMSE
calculate_overall_rmse <- function(predictions) {
  total_sse <- 0
  n_comparisons <- 0
  
  for (site in names(predictions)) {
    site_data <- predictions[[site]]
    differences <- site_data$mean - site_data$observed
    total_sse <- total_sse + sum(differences^2)
    n_comparisons <- n_comparisons + length(differences)
  }
  
  return(sqrt(total_sse / n_comparisons))
}


# Enhanced example usage function with better reporting
example_enhanced_calibration_workflow <- function(patient_data, optimized_params) {
  
  # Enhanced example workflow with comprehensive calibration
  # 
  # Args:
  #   patient_data: Your d4s dataset
  #   optimized_params: Your optimized parameters dataframe
  
  cat("=== Starting Enhanced Calibration Workflow ===\n\n")
  
  # Step 1: Get initial predictions
  cat("Step 1: Getting initial predictions...\n")
  initial_predictions <- predict_stage_distribution(patient_data, optimized_params, n_sim = 1000)
  initial_rmse <- calculate_overall_rmse(initial_predictions)
  cat("Initial RMSE:", round(initial_rmse, 4), "\n\n")
  
  # Step 2: Run enhanced iterative calibration
  cat("Step 2: Running enhanced iterative calibration...\n")
  calibration_results <- iterative_calibration(
    patient_data = patient_data,
    initial_params = optimized_params,
    max_iterations = 8,
    convergence_threshold = 0.015,
    adjustment_weight = 0.25,
    n_sim = 1000,
    adaptive_weight = TRUE
  )
  
  # Step 3: Validate final parameters
  cat("Step 3: Validating final parameters...\n")
  validation_results <- validate_calibrated_parameters(calibration_results$calibrated_params)
  print(validation_results)
  
  # Step 4: Return comprehensive results
  return(list(
    original_predictions = initial_predictions,
    calibration_results = calibration_results,
    # comparison_plot = comparison_plot,
    validation_results = validation_results
  ))
}


# Function to validate all calibrated parameters
validate_calibrated_parameters <- function(calibrated_params) {
  
  validation_summary <- data.frame()
  
  for (i in 1:nrow(calibrated_params)) {
    site <- calibrated_params$site[i]
    params <- as.numeric(calibrated_params[i, c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")])
    
    is_valid <- validate_single_site_parameters(params[1], params[2], params[3], 
                                                params[4], params[5], params[6], params[7])
    
    # Check individual constraints
    bounds_ok <- all(params >= 0) && all(params <= 1)
    sum_constraints_ok <- (params[1] + params[4]) <= 1.001 && 
      (params[2] + params[5]) <= 1.001 && 
      (params[3] + params[6]) <= 1.001 && 
      params[7] <= 1.001
    diagnosis_ordering_ok <- params[4] <= params[5] * 1.1 && 
      params[5] <= params[6] * 1.1 && 
      params[6] <= params[7] * 1.1
    
    validation_summary <- rbind(validation_summary, data.frame(
      site = site,
      overall_valid = is_valid,
      bounds_ok = bounds_ok,
      sum_constraints_ok = sum_constraints_ok,
      diagnosis_ordering_ok = diagnosis_ordering_ok
    ))
  }
  
  return(validation_summary)
}