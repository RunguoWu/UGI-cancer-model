# This is the function file for parameter estimation
# Use different methods of parameter search
# Maximum likelihood estimation for individual patient data

# Search methods -----------------------------------------------------------

# Modified estimate_parameters function using random search
# The mostly used method according to Zhang et al.2025
# DOI: 10.1177/0272989X251353211
estimate_parameters <- function(patient_data,
                                n_sim = 10000, # repetitions for predicting initial distribution
                                n_core = 1,
                                seed_n =1234, 
                                n_cycles = 24,
                                ### Whether consider month in calculating likelihood
                                calib_use_month = TRUE,
                                ###################
                                ### boundries for suggesting starting parameters 
                                dw_list,
                                ####################
                                ### Whether add age, sex and ethnicity in parameters
                                plus_param = FALSE, # TRUE = plus age sex and ethn
                                # Boundries for randomly suggesting parameters
                                suggest_logit_mean = c(-0.16, -0.19, -0.09),
                                suggest_logit_sd = c(0.5, 0.5, 0.5),
                                lower = -1.5,
                                upper = 0.2,
                                sigma = 5,
                                ####################
                                optim_method = "RS", # RS for Random Search; NM for Nelder-Mead; LBB for L-BFGS-B
                                n_param_try # Number of random search or initial param for other methods
                                ) {
  
  set.seed(seed_n)

  # Objective function (negative log-likelihood) with progress tracking
  objective <- function(params) {
    
    # Validate parameters first
    if (!validate_parameters(params)) {
      return(1e8)  # Return large but finite penalty
    }
    
    # Try to calculate likelihood with error handling
    likelihood <- tryCatch({
      calculate_cohort_likelihood(params, patient_data, 
                                  use_month = calib_use_month,
                                  plus_param = plus_param,
                                  n_cycles = n_cycles,
                                  n_sim = n_sim,
                                  sigma = sigma
                                  )
    }, error = function(e) {
      cat("Parameters:", params, "\n")
      cat("Error:", e$message, "\n")
      return(-Inf)
    })
    
    # Handle infinite likelihood
    if (is.infinite(likelihood)) {
      return(1e9)  # Return large but finite penalty
    }
    
    # Handle NaN or NA
    if (is.na(likelihood)) {
      return(1e10)  # Return large but finite penalty
    }
    
    return(-likelihood)
  }
  
  method_name <- case_when(optim_method == "RS" ~ "Random Search",
                           optim_method == "NM" ~ "Nelder-Mead",
                           optim_method == "LBB" ~ "L-BFGS-B", 
                           TRUE ~ "ERROR: incorrect method"
                           )
  
  if (method_name == "ERROR: incorrect method") {
    warning(method_name)
    return(NULL)
  }
   
  cat("Starting optimization for 4-stage model using", method_name, "...\n")
  cat("Number of initial samples:", n_param_try, "\n")
  
  if (optim_method == "LBB"){
    bounds <- get_parameter_bounds(plus_param)
    lower_params <- bounds$lower
    upper_params <- bounds$upper
    
    cat("Parameter bounds: [", lower_params, "] to [", upper_params, "]\n")
  }
  
  stage_prop <- prop.table(table(patient_data$diagnosed_stage))

  # Pre-sampling of parameters
  random_params_list <- lapply(1:n_param_try, function(i) {
    suggest_starting_parameters(stage_prop,
                                dw_list,
                                plus_param = plus_param,
                                suggest_logit_mean = suggest_logit_mean,
                                suggest_logit_sd = suggest_logit_sd,
                                lower = lower,
                                upper = upper)
  })
  
  random_params_list <- do.call(rbind, random_params_list)
  
  if (optim_method == "RS"){
    # Random search optimization
    # Initialize with initial parameters
    initial_params <- random_params_list[1, ]
    best_params <- initial_params
    best_value <- objective(initial_params)
    
    cat("Initial objective value:", best_value, "\n\n")
  }
  
  start_time <- Sys.time()
  
  cl <- makePSOCKcluster(n_core)
  clusterExport(cl, c(# "random_params_list", "objective", "best_params", "best_value", 
                      "validate_parameters", "subset_data",
                      "create_transition_matrix", "simulate_markov", 
                      "initial_stage_distribution", "save_intermediate_results", 
                      "calculate_cohort_likelihood",
                      "tp_to_rate"
                      ))
  clusterEvalQ(cl, library(tidyverse))
  registerDoSNOW(cl)

  if (optim_method == "RS") {
    # Generate random parameter sets
    # for (i in 1:n_param_try)  {
    best_list <- foreach (i = 2:n_param_try) %dopar% {
      # Generate random parameters within bounds
      random_params <- random_params_list[i, ]
      
      # Evaluate objective function
      current_value <- objective(random_params)
      
      # Update best if current is better (lower negative log-likelihood)
      if (current_value < best_value) {
        best_params <- random_params
        best_value <- current_value
      }
      
      # Progress update every 1000 iterations
      if (i %% 1000 == 0) {
        cat("Completed", i, "of", n_param_try, "samples\n")
      }
      
      rt <- list(value = best_value, par = best_params, 
                 value_cur = current_value, par_cur = random_params)
      return(rt)
    }
    stopCluster(cl)
  }
  
  if (optim_method == "NM") {
    best_list <- foreach (i = 1:n_param_try) %dopar% {
    
      initial_params <- random_params_list[i, ]

      rt <- optim(
        par = initial_params,
        fn = objective,
        method = "Nelder-Mead",
        control = list(
          trace = 1,
          maxit = 6000,
          reltol = 1e-8)
      )
      return(rt)
    }
    stopCluster(cl)
  }
  
  if (optim_method == "LBB") {
    
    best_list <- foreach (i = 1:n_param_try) %dopar% {
      
      initial_params <- random_params_list[i, ]
      
      rt <- optim(
        par = initial_params,
        fn = objective,
        method = "L-BFGS-B",
        lower = lower_params,
        upper = upper_params,
        control = list(
          trace = 1,
          maxit = 6000,
          factr = 1e12,
          pgtol = 1e-8,
          ndeps = rep(1e-8, length(initial_params))
        )
      )
      return(rt)
    }
    stopCluster(cl)
  }
  
  end_time <- Sys.time()
  
  if (optim_method == "NM") {
    # As Nelder-Mead method does not have bounds
    # use validate function to screen the parameters first.
    best_list <- best_list[sapply(best_list, function(x) validate_parameters(x$par))]
  }
  
  best_value_list <- sapply(best_list, "[[", "value")
  best_value <- best_list[[which.min(best_value_list)]][["value"]]
  best_params <- best_list[[which.min(best_value_list)]][["par"]]

  if (optim_method == "RS") {
    params_record <- lapply(best_list, function(x) x[c("value_cur", "par_cur")])
  } else params_record <- NULL
  
  if (optim_method == "NM") {
    best_result <- best_list[[which.min(best_value_list)]]
    convergence_code  <-  best_result$convergence
    convergence_message  <-  switch(as.character(convergence_code),
                                 "0" = "Successful convergence",
                                 "1" = "Maximum iterations reached",
                                 "10" = "Degenerate simplex",
                                 paste("Unknown convergence code:", best_result$convergence))
  }
  
  if (optim_method == "RS") {
    convergence_code <- 0
    convergence_message <- paste("Random search completed with", n_param_try, "samples")
  }
  
  end_time <- Sys.time()
  
  # Validate final results
  if (!validate_parameters(best_params)) {
    warning("Final parameters violate constraints. Results may not be reliable.")
  }
  
  parameter_names <- c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
  if (plus_param) parameter_names <- c(parameter_names, "beta_age70plus", "beta_female", "beta_nonwhite")
  
  # Save final results
  final_result <- list(
    value_convergence_record = best_value_list,
    params_record = params_record,
    computation_time = end_time - start_time,
    final_params = best_params,
    final_likelihood = -best_value,
    convergence_code = convergence_code,
    convergence_message = convergence_message,
    parameter_names = parameter_names,
    n_param_try = n_param_try
  )
  
  cat("\nOptimization complete!\n")
  cat("Method:", method_name, "\n")
  cat("Convergence:", final_result$convergence_message, "\n")
  cat("Total samples evaluated:", n_param_try, "\n")
  cat("Total time:", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")
  
  cat("Final parameters:\n")
  for (i in 1:length(best_params)) {
    cat(sprintf("  %s: %.4f\n", final_result$parameter_names[i], best_params[i]))
  }
  cat("Final likelihood:", round(-best_value, 4), "\n")
  
  return(final_result)
}


# Recalibration -----------------------------------------------------------

estimate_parameters_recali <- function(patient_data,
                                    optimized_params,
                                    filename_base = "mle_results_recali",
                                    save_final = FALSE, 
                                    seed_n =1234
) {
  
  set.seed(seed_n)
  
  # Objective function (negative log-likelihood) with progress tracking
  objective <- function(params) {
    
    # Validate parameters first
    if (!validate_parameters(params)) {
      return(1e7)  # Return large but finite penalty
    }
    
    if (any(params < optimized_params * 0.2) || any(params > optimized_params * 5)) return(1e10)
    
    if (params[2] < params[1]) return(1e10)  
    if (params[3] < params[2]) return(1e10)  
    
    # 2. Detection rates should increase with stage
    if (params[5] < params[4]) return(1e10)  
    if (params[6] < params[5]) return(1e10)  
    if (params[7] < params[6]) return(1e10) 
    
    # Try to calculate likelihood with error handling
    likelihood <- tryCatch({
      calculate_cohort_likelihood(params, patient_data, use_month = FALSE)
    }, error = function(e) {
      cat("Parameters:", params, "\n")
      cat("Error:", e$message, "\n")
      return(NA)
    })
    
    # Handle infinite likelihood
    if (is.infinite(likelihood)) {
      return(1e8)  # Return large but finite penalty
    }
    
    # Handle NaN or NA
    if (is.na(likelihood)) {
      return(1e9)  # Return large but finite penalty
    }
    
    return(-likelihood)
  }
  
  site_val <- unique(patient_data$site)
  
  cat("Starting optimization for 4-stage model...\n", site_val)
  
  start_time <- Sys.time()
  
  # Optimization with constraints
  # Run optimParallel
  # rt <- optim(
  #   par = optimized_params,
  #   fn = objective,
  #   method = "L-BFGS-B",
  #   lower = optimized_params*0.2,
  #   upper = optimized_params*5,
  #   control = list(
  #     trace = 1,
  #     maxit = 6000,
  #     factr = 1e15,
  #     parscale = abs(optimized_params) + 0.1
  #   )
  # )
  # # 
  # rt <- optim(optimized_params, objective, method = "L-BFGS-B",
  #       lower = optimized_params * 0.5,
  #       upper = optimized_params * 2,
  #       control = list(factr = 1e12, lmm = 15))
  
  rt0 <- optim(optimized_params, objective, method = "Nelder-Mead",
               control = list(
                 reltol = 1e-8,
                 maxit = 5000
               ))
  

  rt <- nlminb(
    start = rt0$par,  # Use Nelder-Mead result
    objective = objective,
    lower = rt0$par / 5,
    upper = rt0$par * 5,
    control = list(trace = 1, iter.max = 1000)
  )
  
  end_time <- Sys.time()
  
  best_params <- rt$par
  
  # Validate final results
  if (!validate_parameters(best_params)) {
    warning("Final parameters violate constraints. Results may not be reliable.")
  }
  
  # Save final results
  final_result <- list(
    mle_result = rt,
    computation_time = end_time - start_time,
    final_params = best_params,
    final_likelihood = -rt$objective,
    convergence_code = rt$convergence,
    convergence_message = switch(as.character(rt$convergence),
                                 "0" = "Successful convergence",
                                 "1" = "Maximum iterations reached",
                                 "51" = "Warning from L-BFGS-B",
                                 "52" = "Error from L-BFGS-B",
                                 paste("Unknown convergence code:", rt$convergence)),
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
  
  return(final_result)
}







