# Calibration Functions for Cancer Natural History Model Parameters

library(tidyverse)

# Enhanced calibration adjustment factors calculation
#######################################################
### This is the key function for the re-calibration ###
#######################################################
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
      progression_adjustments[1] <- progression_adjustments[1] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.7)
      progression_adjustments[2] <- progression_adjustments[2] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.5)
      progression_adjustments[3] <- progression_adjustments[3] * (1 + target_weight * (obs_pred_ratios[4] - 1) * 0.3)
    }
    
    # Diagnosis rate adjustments (direct from obs/pred ratios)
    # Use as if under-predicted happens to stage 1-3
    # tp4 does not affect stage distribution at all
    diagnosis_adjustments[1:3] <- 1 + target_weight* 2 * (obs_pred_ratios[1:3] - 1) # calibrate diagnosis rate faster
    
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
    
    site_row <- which(calibrated_params$site == site)
    diagnosis_factors <- calibration_factors[[site]]$diagnosis
    progression_factors <- calibration_factors[[site]]$progression
    
    # Store original values for comparison
    original_values <- as.numeric(calibrated_params[site_row, c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")])
    
    # # Apply factors to progression rate parameters (tp12, tp23, tp34)
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
    if (progression_rates[stage] > progression_rates[stage + 1]) {
      # Only adjust if the violation is significant
      avg_rate <- (progression_rates[stage] + progression_rates[stage + 1]) / 2
      progression_rates[stage] <- avg_rate * 0.95
      progression_rates[stage + 1] <- avg_rate * 1.05
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
                                  max_iterations = 15, 
                                  convergence_threshold = 0.015,
                                  adjustment_weight = 0.25,
                                  n_sim = 100,
                                  months = 24,
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
  baseline_predictions <- predict_stage_distribution(patient_data, 
                                                     current_params, 
                                                     n_sim = n_sim, months = months)
  baseline_rmse <- calculate_overall_rmse(baseline_predictions)
  
  repeat {
    iteration <- iteration + 1
    cat("=== Calibration Iteration", iteration, "===\n")
    
    # Predict with current parameters
    predictions <- predict_stage_distribution(patient_data, 
                                              current_params, 
                                              n_sim = n_sim, months = months)
    
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
  final_predictions <- predict_stage_distribution(patient_data, current_params, 
                                                  n_sim = n_sim, months = months)
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
  initial_predictions <- predict_stage_distribution(patient_data, 
                                                    optimized_params, 
                                                    n_sim = 100, months = 24)
  initial_rmse <- calculate_overall_rmse(initial_predictions)
  cat("Initial RMSE:", round(initial_rmse, 4), "\n\n")
  
  # Step 2: Run enhanced iterative calibration
  cat("Step 2: Running enhanced iterative calibration...\n")
  calibration_results <- iterative_calibration(
    patient_data = patient_data,
    initial_params = optimized_params,
    max_iterations = 10,
    convergence_threshold = 0.015,
    adjustment_weight = 0.25,
    n_sim = 100,
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