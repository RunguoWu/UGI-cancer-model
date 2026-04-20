# Functions of optimisation of natural history parameters - 2-stage version

library(tidyverse)
library(expm)  # for matrix exponentiation

# Function to create transition matrix
# state:
# 1: undiagnosed early
# 2: undiagnosed late
# 3: diagnosed early
# 4: diagnosed late

create_transition_matrix <- function(tp1, tp2, tp3) {
  # tp1: undiagnosed early -> undiagnosed late
  # tp2: undiagnosed early -> diagnosed early  
  # tp3: undiagnosed late -> diagnosed late
  
  # Probability of staying in undiagnosed early
  stay_early <- 1 - tp1 - tp2
  
  # Probability of staying in undiagnosed late
  stay_late <- 1 - tp3
  
  # Create transition matrix
  P <- matrix(c(
    stay_early, tp1,       tp2, 0,    # From undiagnosed early
    0,          stay_late, 0,   tp3,  # From undiagnosed late
    0,          0,         1,   0,    # From diagnosed early (absorbing)
    0,          0,         0,   1     # From diagnosed late (absorbing)
  ), nrow = 4, byrow = TRUE)
  
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
# initial_state = 1: undiagnosed early stage, column 1 in matrix
# initial_state = 2: undiagnosed late stage, column 2 in matrix
simulate_markov <- function(P, initial_state, n_cycles = 24) {
  # Initialize state distribution
  state_dist <- rep(0, 4)
  state_dist[initial_state] <- 1
  
  # Store results for each cycle
  results <- matrix(0, nrow = n_cycles + 1, ncol = 4)
  results[1, ] <- state_dist
  
  # Simulate each cycle
  for (cycle in 1:n_cycles) {
    state_dist <- state_dist %*% P
    results[cycle + 1, ] <- state_dist
  }
  
  return(results)
}

# Function to calculate probability of starting in early stage
# I do not optimise the parameter directly, as it is likely to go to the lowest bound
# to maximise late stage start, which is easy to get better likelihood
# as it has only two stages to guess: undiagnosed late and diagnosed late
# Instead, I use tp1 and tp3 to guess the probability 
# this function only applies to diagnosed at late stage
calculate_early_start_probability <- function(month_diagnosis, tp1, tp3, 
                                              n_sim = 10000, seed_num=1234) {
  set.seed(seed_num)
  time_early_late <- rexp(n_sim, rate = tp1)
  time_late_diag <- rexp(n_sim, rate = tp3)
  
  # Vectorized conditions
  # time from index to dignosis > total time from late to diagnosis, suggesting 
  # there is remaining time for transition from early to late before diagnosis
  # conditioning on remaining time is not longer than the total time from early to late
  condition1 <- (time_late_diag < month_diagnosis) & 
    ((month_diagnosis - time_late_diag) < time_early_late)
  
  # time from index to dignosis < total time from late to diagnosis, suggesting
  # transition from early to late happened before index
  condition2 <- time_late_diag > month_diagnosis
  
  early_start_count <- sum(condition1)
  late_start_count <- sum(condition2)
  
  # Handle edge case where no valid transitions occur
  if (early_start_count + late_start_count == 0) {
    return(0.5)  # Return neutral probability
  }
  
  prob_early_start <- early_start_count / (early_start_count + late_start_count)
  
  return(prob_early_start)
}

# Function to calculate total log-likelihood for individual patient cohort
calculate_cohort_likelihood <- function(params, patient_data) {
  # patient_data should be a data frame with columns:
  # - month: month of diagnosis (0-24)
  # - diagnosed_stage: 3 for early, 4 for late
  
  # Input validation
  if (length(params) != 3) return(-Inf)
  
  if (any(is.na(params)) || any(is.null(params))) return(-Inf)
  
  # Extract parameters
  tp1 <- params[1]  # undiagnosed early -> undiagnosed late
  tp2 <- params[2]  # undiagnosed early -> diagnosed early
  tp3 <- params[3]  # undiagnosed late -> diagnosed late
  
  # Parameter constraints
  if (tp1 < 0 || tp1 > 1 || tp2 < 0 || tp2 > 1 || tp3 < 0 || tp3 > 1 ||
      (tp1 + tp2) > 1) {
    return(-Inf)
  }
  
  # Create transition matrix
  P <- create_transition_matrix(tp1, tp2, tp3)
  
  # Pre-compute simulations to avoid repeated calculations
  sim_early <- simulate_markov(P, initial_state = 1, n_cycles = 24)
  sim_late <- simulate_markov(P, initial_state = 2, n_cycles = 24)
  
  patient_summary <- patient_data %>%
    group_by(month, diagnosed_stage) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate total log-likelihood
  total_log_likelihood <- 0
  
  for (i in 1:nrow(patient_summary)) {
    month <- patient_summary$month[i]
    diagnosed_stage <- patient_summary$diagnosed_stage[i]
    count <- patient_summary$count[i]
    
    if (month < 0 || month > 24) next
    
    if (diagnosed_stage == 3) {
      # if diagnosed at early stage, we assume that the patient was at early stage at index date
      # Early diagnosis case
      if (month + 1 <= nrow(sim_early)) {
        if (month == 0) {
          # At month 0, use cumulative probability of early diagnosis
          prob <- sim_early[month + 1, 3]
        } else {
          # For month > 0, we want P(undiag early at month-1 AND diag early at month)
          # Method 1: Using transition probability
          prob_undiag_early_prev <- sim_early[month, 1]
          prob <- prob_undiag_early_prev * tp2
        }
        
        prob <- max(prob, 1e-10)
        total_log_likelihood <- total_log_likelihood + count * log(prob)
      }
      
    } else if (diagnosed_stage == 4) {
      # if diagnosed at early stage, 
      # we consider both possibilities that the patient was at early or late stage on the index date
      # Late diagnosis - use probabilistic weighting
      prob_early_start <- calculate_early_start_probability(month, tp1, tp3)
      
      # Check if early start probability is valid
      if (is.na(prob_early_start) || is.infinite(prob_early_start)) {
        return(-Inf)
      }
      
      # Calculate likelihood for both scenarios
      if (month + 1 <= nrow(sim_early)) {
        # Likelihood if started in early stage
        if (month == 0) {
          prob_early_scenario <- sim_early[month + 1, 4]
        } else {
          prob_undiag_late_prev <- sim_early[month, 2]
          prob_early_scenario <- prob_undiag_late_prev * tp3
        }
        
        # Likelihood if started in late stage
        if (month == 0) {
          prob_late_scenario <- sim_late[month + 1, 4]
        } else {
          prob_undiag_late_prev <- sim_late[month, 2]
          prob_late_scenario <- prob_undiag_late_prev * tp3
        }
        
        # Weighted mixture
        prob <- prob_early_start * prob_early_scenario + 
          (1 - prob_early_start) * prob_late_scenario
        
        prob <- max(prob, 1e-10)
        total_log_likelihood <- total_log_likelihood + count * log(prob)
      }
    }
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
  
  saveRDS(result, file = file.path(output,"optim", filename))
  cat("Results saved to:", filename, "\n")
  return(filename)
}

# Robust parameter validation function
validate_parameters <- function(params) {
  if (length(params) != 3) return(FALSE)
  if (any(is.na(params)) || any(is.null(params))) return(FALSE)
  
  tp1 <- params[1]
  tp2 <- params[2]
  tp3 <- params[3]
  
  # Check bounds
  if (any(params < 0) || any(params > 1)) return(FALSE)
  
  # Check constraint: tp1 + tp2 <= 1
  if (tp1 + tp2 > 1) return(FALSE)
  
  # Check constraint: tp2 < tp3
  if (tp2 >= tp3) return(FALSE)
  
  return(TRUE)
}


# Maximum likelihood estimation for individual patient data
estimate_parameters <- function(patient_data, initial_params, 
                                lower_params, upper_params,
                                save_interval = NA,
                                filename_base = "mle_results",
                                save_final = TRUE
) {
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
    # save_interval = NA: do not save intermediate results
    if (!is.na(save_interval) & iteration_count %% save_interval == 0) {
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
  
  cat("Starting optimization...\n")
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
        maxit = 2000, 
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
      control = list(trace = 1, maxit = 1000)
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
                                 paste("Unknown convergence code:", result$convergence))
  )
  
  final_filename <- NULL
  if (save_final){
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
  cat("Final parameters:", round(result$par, 4), "\n")
  cat("Final likelihood:", round(-result$value, 4), "\n")
  
  if (!is.null(final_filename)) {
    cat("Final results saved to:", final_filename, "\n")
  }
  
  return(final_result)
}

# Helper function to suggest good starting parameters
suggest_starting_parameters <- function(patient_data) {
  # Basic analysis of the data
  early_count <- sum(patient_data$diagnosed_stage == 3)
  late_count <- sum(patient_data$diagnosed_stage == 4)
  total_count <- nrow(patient_data)
  
  # Calculate proportion diagnosed early vs late
  early_prop <- early_count / total_count
  late_prop <- late_count / total_count
  
  # Suggest parameters based on data characteristics
  suggested_tp1 <- 0.1  # Conservative progression rate
  suggested_tp2 <- max(0.01, early_prop * 0.2)  # Based on early diagnosis rate
  suggested_tp3 <- max(suggested_tp2 * 1.5, late_prop * 0.3)  # Higher than tp2
  
  cat("Data summary:\n")
  cat("Total patients:", total_count, "\n")
  cat("Early stage:", early_count, "(", round(early_prop * 100, 1), "%)\n")
  cat("Late stage:", late_count, "(", round(late_prop * 100, 1), "%)\n")
  cat("\nSuggested starting parameters:\n")
  cat("tp1 (early->late progression):", round(suggested_tp1, 4), "\n")
  cat("tp2 (early->early diagnosis):", round(suggested_tp2, 4), "\n")
  cat("tp3 (late->late diagnosis):", round(suggested_tp3, 4), "\n")
  
  return(c(suggested_tp1, suggested_tp2, suggested_tp3))
}
