# Validation 2, microsimulation -------------------------------------------

# Function to estimate starting stage distribution for each cancer site
estimate_starting_distribution <- function(d4s, params, plus_param=FALSE, average_dist = FALSE) {
  # Group by cancer site and calculate stage distribution for each
  sites <- unique(d4s$site)
  starting_distributions <- list()
  
  for (st in sites) {
    site_data <- d4s[d4s$site == st, ]
    
    optimized_params <- params[params$site==st, grepl("tp|beta", colnames(params))]
    
    tp12 <- as.numeric(optimized_params["tp12"])  
    tp23 <- as.numeric(optimized_params["tp23"]) 
    tp34 <- as.numeric(optimized_params["tp34"]) 
    tp1 <- as.numeric(optimized_params["tp1"])   
    tp2 <- as.numeric(optimized_params["tp2"]) 
    tp3 <- as.numeric(optimized_params["tp3"])   
    tp4 <- as.numeric(optimized_params["tp4"])   
    
    if (plus_param){
      beta_age70plus <- as.numeric(optimized_params["beta_age70plus"])
      beta_female <- as.numeric(optimized_params["beta_female"])
      beta_nonwhite <- as.numeric(optimized_params["beta_nonwhite"])
    }
    
    # Calculate weighted average of starting distributions across all patients
    stage_dist_matrix <- matrix(0, nrow = nrow(site_data), ncol = 4)
    
    for (i in 1:nrow(site_data)) {
      
      if (plus_param){
        female <- site_data$female[i]
        nonwhite <- site_data$nonwhite[i]
        age70plus <- site_data$age70plus[i]
        
        #   multiplier <- exp(beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite)
        # } else {multiplier <- 1}
        
        multiplier <- beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite
      } else {multiplier <- 0}
        
      stage_dist_matrix[i, ] <- initial_stage_distribution(
        month_diagnosis = site_data$month[i],
        diagnosed_stage = site_data$diagnosed_stage[i],
        tp12 = tp12, tp23 = tp23, tp34 = tp34,
        tp1 = tp1, tp2 = tp2, tp3 = tp3, tp4 = tp4, 
        multiplier = multiplier
      )
    }
    
    stage_dist_df <- data.frame(
      e_patid = as.character(site_data$e_patid),
      stage_dist_matrix,
      stringsAsFactors = FALSE
    )
    
    # Average across patients for this cancer site
    starting_distributions[[st]] <- if(average_dist) colMeans(stage_dist_matrix) else stage_dist_df
  }
  
  return(starting_distributions)
}

# Microsimulation function for single patient
simulate_single_patient <- function(P, starting_stage_probs, months=60) {
  # Sample starting stage based on probabilities
  start_stage <- sample(1:4, 1, prob = starting_stage_probs)
  
  current_state <- start_stage
  
  for (cycle in 1:months) {

    # Sample next state based on transition probabilities
    next_state <- sample(1:8, 1, prob = P[current_state, ])
    current_state <- next_state
    
    # If transitioned to diagnosed state, return diagnosis stage
    if (current_state > 4) {
      return(current_state - 4)
    }
    if (current_state == 4) {
      return(current_state)
    }
  }
  
  # If still undiagnosed after 24 cycles, return current stage
  # return(current_state)
  
  # If still undiagnosed after months cycles, return NA
  return(NA)
}

# Main prediction function
# Cohort level
# The result is same as patient level, but
# this is faster, so use it for calibration
predict_stage_distribution <- function(patient_data, params, n_sim = 500, 
                                       months=60, plus_param = FALSE) {
  
  # Estimate starting distributions by cancer site
  starting_distributions <- estimate_starting_distribution(patient_data, params, plus_param=plus_param)
  
  # Initialize results storage
  sites <- names(starting_distributions)
  predicted_distributions <- list()
  
  if (plus_param){
    # age70plus == 1 & female == 1 & nonwhite == 1, "o7fn"
    # age70plus == 0 & female == 0 & nonwhite == 0, "u7mw" 
    addon_names <- c("o7_fe_nw", "o7_fe_wh", "o7_ma_nw", "o7_ma_wh", "u7_fe_nw", "u7_fe_wh", "u7_ma_nw", "u7_ma_wh")
    addon_names_bi <- c("o7", "u7", "fe", "ma", "nw", "wh")
  }
  
  for (st in sites) {
    site_results <- matrix(0, nrow = n_sim, ncol = 4)
    diagnosis_rates <- numeric(n_sim)  # Track proportion diagnosed
    site_data <- patient_data[patient_data$site == st, ]
    n_patients <- nrow(site_data)
    
    # Create transition matrix
    optimized_params <- params[params$site==st, grepl("tp|beta", colnames(params))]
    
    tp12 <- as.numeric(optimized_params["tp12"])  
    tp23 <- as.numeric(optimized_params["tp23"]) 
    tp34 <- as.numeric(optimized_params["tp34"]) 
    tp1 <- as.numeric(optimized_params["tp1"])   
    tp2 <- as.numeric(optimized_params["tp2"]) 
    tp3 <- as.numeric(optimized_params["tp3"])   
    tp4 <- as.numeric(optimized_params["tp4"])   
    
    if (plus_param){
      beta_age70plus <- as.numeric(optimized_params["beta_age70plus"])
      beta_female <- as.numeric(optimized_params["beta_female"])
      beta_nonwhite <- as.numeric(optimized_params["beta_nonwhite"])
      
      site_results_sub <- list()
      for (name in addon_names) {
        site_results_sub[[paste0("site_results_", name)]] <- matrix(0, nrow = n_sim, ncol = 4)
      }
      
      for (name in addon_names_bi) {
        site_results_sub[[paste0("site_results_", name)]] <- matrix(0, nrow = n_sim, ncol = 4)
      }
    }
    
    P <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
    
    cat("Simulating", st, "cancer (", n_patients, "patients)...\n")
    
    for (sim in 1:n_sim) {
      
      counts_list <- list(
        stage_counts = rep(0, 4),
        n_diagnosed = 0
      )
      
      # count subgroups
      if (plus_param){
        for (name in addon_names) {
          counts_list[[paste0("stage_counts_", name)]] <- rep(0, 4)
          counts_list[[paste0("n_diagnosed_", name)]] <- 0
        }
      }

      # Simulate each patient in this cancer site
      for (i in 1:n_patients) {
        
        if (plus_param){
          female <- site_data$female[i]
          nonwhite <- site_data$nonwhite[i]
          age70plus <- site_data$age70plus[i]
          
          # multiplier <- exp(beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite)
          
          multiplier <- beta_age70plus*age70plus + beta_female*female + beta_nonwhite*nonwhite
          
          PP <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4, multiplier = multiplier)
          
          # pick up one of these addon names
          addon_name <- case_when(
            age70plus == 1 & female == 1 & nonwhite == 1 ~ "o7_fe_nw",  # 1,1,1
            age70plus == 1 & female == 1 & nonwhite == 0 ~ "o7_fe_wh",  # 1,1,0
            age70plus == 1 & female == 0 & nonwhite == 1 ~ "o7_ma_nw",  # 1,0,1
            age70plus == 1 & female == 0 & nonwhite == 0 ~ "o7_ma_wh",  # 1,0,0
            age70plus == 0 & female == 1 & nonwhite == 1 ~ "u7_fe_nw",  # 0,1,1
            age70plus == 0 & female == 1 & nonwhite == 0 ~ "u7_fe_wh",  # 0,1,0
            age70plus == 0 & female == 0 & nonwhite == 1 ~ "u7_ma_nw",  # 0,0,1
            age70plus == 0 & female == 0 & nonwhite == 0 ~ "u7_ma_wh",  # 0,0,0
            TRUE ~ NA_character_  # For any missing values
          )
        } else {
          PP <- P
        }
        
        diagnosed_stage <- simulate_single_patient(PP, starting_distributions[[st]][i, -1], months=months) # [i, -1] remove the id, only keep 4 stages
        
        if (!is.na(diagnosed_stage)) {
          counts_list[["stage_counts"]][diagnosed_stage] <- counts_list[["stage_counts"]][diagnosed_stage] + 1
          counts_list[["n_diagnosed"]] <- counts_list[["n_diagnosed"]] + 1
          
          if (plus_param){
            counts_list[[paste0("stage_counts_", addon_name)]][diagnosed_stage] <- 
              counts_list[[paste0("stage_counts_", addon_name)]][diagnosed_stage] + 1
            
            counts_list[[paste0("n_diagnosed_", addon_name)]] <- 
              counts_list[[paste0("n_diagnosed_", addon_name)]] + 1
          }
        }
      }
      
      # Store proportions for this simulation (among diagnosed patients only)
      if (counts_list[["n_diagnosed"]] > 0) {
        site_results[sim, ] <- counts_list[["stage_counts"]] / counts_list[["n_diagnosed"]] # n_patients
      } else {
        site_results[sim, ] <- rep(NA, 4)
      }
    
      if (plus_param){
        
        for (name in addon_names) {
          
          if (counts_list[[paste0("n_diagnosed_", name)]] > 0){
            site_results_sub[[paste0("site_results_", name)]][sim, ] <- 
              counts_list[[paste0("stage_counts_", name)]] / counts_list[[paste0("n_diagnosed_", name)]]
          } else {
            site_results_sub[[paste0("site_results_", name)]][sim, ] <- rep(NA, 4)
          }
        }
        
        # by one characteristic
        for (name in addon_names_bi){
          
          sub_list <- counts_list[grepl(name, names(counts_list)) & grepl("stage_counts", names(counts_list))]
          sub_stage_counts <- colSums(do.call(rbind, sub_list))
          
          sub_list2 <- counts_list[grepl(name, names(counts_list)) & grepl("n_diagnosed", names(counts_list))]
          sub_n_diagnosed <- sum(do.call(c, sub_list2))
          
          if (sub_n_diagnosed > 0){
            site_results_sub[[paste0("site_results_", name)]][sim, ] <- 
              sub_stage_counts / sub_n_diagnosed
          } else {
            site_results_sub[[paste0("site_results_", name)]][sim, ] <- rep(NA, 4)
          }
        }
      }
      
      diagnosis_rates[sim] <- counts_list[["n_diagnosed"]] / n_patients
    }
      
    # Calculate mean and confidence intervals across simulations
    predicted_distributions[[st]] <- list(
      mean = colMeans(site_results, na.rm = TRUE),
      ci_lower = apply(site_results, 2, quantile, 0.025, na.rm = TRUE),
      ci_upper = apply(site_results, 2, quantile, 0.975, na.rm = TRUE),
      observed = as.numeric(table(factor(site_data$diagnosed_stage, levels = 1:4)) / nrow(site_data)),
      observed_lower = obs_ci(site_data)[["obs_lower"]],
      observed_upper = obs_ci(site_data)[["obs_upper"]],
      diagnosis_rate = mean(diagnosis_rates),
      diagnosis_rate_ci = quantile(diagnosis_rates, c(0.025, 0.975))
    )
    
    if (plus_param){
      
      conditions_map <- list(
        o7_fe_nw = quote(age70plus == 1 & female == 1 & nonwhite == 1),
        o7_fe_wh = quote(age70plus == 1 & female == 1 & nonwhite == 0),
        o7_ma_nw = quote(age70plus == 1 & female == 0 & nonwhite == 1),
        o7_ma_wh = quote(age70plus == 1 & female == 0 & nonwhite == 0),
        u7_fe_nw = quote(age70plus == 0 & female == 1 & nonwhite == 1),
        u7_fe_wh = quote(age70plus == 0 & female == 1 & nonwhite == 0),
        u7_ma_nw = quote(age70plus == 0 & female == 0 & nonwhite == 1),
        u7_ma_wh = quote(age70plus == 0 & female == 0 & nonwhite == 0)
      )
      
      # Calculate observed distributions for all subgroups
      # all cells
      observed_list <- list()
      for (name in addon_names) {
        subgroup_data <- subset(site_data, eval(conditions_map[[name]]))
        observed_list[[name]] <- as.numeric(table(factor(subgroup_data$diagnosed_stage, levels = 1:4)) / 
                                              nrow(subgroup_data))
      }
      
      for (name in addon_names) {
        
        predicted_distributions[[st]][["sub_group"]][[name]][["mean"]] <-  
          colMeans(site_results_sub[[paste0("site_results_", name)]], na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["ci_lower"]] <-  
          apply(site_results_sub[[paste0("site_results_", name)]], 2, quantile, 0.025, na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["ci_upper"]] <-  
          apply(site_results_sub[[paste0("site_results_", name)]], 2, quantile, 0.975, na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["observed"]] <- observed_list[[name]]
      }
      
      # by single characteristic
      conditions_map2 <- list(
        o7 = quote(age70plus == 1),
        fe = quote(female == 1),
        nw = quote(nonwhite == 1),
        u7 = quote(age70plus == 0),
        ma = quote(female == 0),
        wh = quote(nonwhite == 0)
      )
      
      # Calculate observed distributions for all subgroups
      observed_list2 <- list()
      for (name in addon_names_bi) {
        subgroup_data <- subset(site_data, eval(conditions_map2[[name]]))
        observed_list2[[name]] <- as.numeric(table(factor(subgroup_data$diagnosed_stage, levels = 1:4)) / 
                                              nrow(subgroup_data))
      }
      
      for (name in addon_names_bi) {
        
        predicted_distributions[[st]][["sub_group"]][[name]][["mean"]] <-  
          colMeans(site_results_sub[[paste0("site_results_", name)]], na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["ci_lower"]] <-  
          apply(site_results_sub[[paste0("site_results_", name)]], 2, quantile, 0.025, na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["ci_upper"]] <-  
          apply(site_results_sub[[paste0("site_results_", name)]], 2, quantile, 0.975, na.rm = TRUE)
        
        predicted_distributions[[st]][["sub_group"]][[name]][["observed"]] <- observed_list2[[name]]
      }
    }
  }
  
  return(predicted_distributions)
}

# Patient level 
# allow treatment effect at the individual level
predict_stage_distribution2 <- function(patient_data, params, n_sim = 500, months=24,
                                        tx, 
                                        tx_prob=1, # proportion receiving tx
                                        use_avg_start_dist = FALSE, 
                                        avg_start_dist = NULL
                                        ) {
  
  # Estimate starting distributions by cancer site
  starting_distributions <- estimate_starting_distribution(patient_data, params)
  
  # Initialize results storage
  sites <- names(starting_distributions)
  predicted_distributions <- list()
  
  for (st in sites) {
    site_data <- patient_data[patient_data$site == st, ]
    n_patients <- nrow(site_data)
    
    # Store results at the simulation level (like original code)
    simulation_results <- matrix(0, nrow = n_sim, ncol = 4)
    simulation_diagnosis_rates <- numeric(n_sim)
    
    # Store patient-level info for additional analysis
    patient_stage_counts <- matrix(0, nrow = n_patients, ncol = 4)
    patient_diagnosed_counts <- numeric(n_patients)
    
    # Create transition matrix
    optimized_params <- as.numeric(params[params$site==st, grepl("tp", colnames(params))])
    P <- create_transition_matrix(
      optimized_params[1], optimized_params[2], optimized_params[3],
      optimized_params[4], optimized_params[5], optimized_params[6], optimized_params[7]
    )
    
    tx_st <- tx[st]
    
    if (tx_st != 1){ # >1 faster detection; <1 slower detection
      P_tx <- create_transition_matrix(
        optimized_params[1], optimized_params[2], optimized_params[3],
        optimized_params[4]*tx_st, optimized_params[5]*tx_st, optimized_params[6]*tx_st, 
        optimized_params[7]
      )
      
      # which patients would receive tx
      tx_patient <- sample(c(TRUE, FALSE), size = n_patients, replace = TRUE, prob = c(tx_prob, 1-tx_prob))
    }
    
    cat("Simulating", st, "cancer (", n_patients, "patients)...\n")
    
    # Loop through each patient first
    for (i in 1:n_patients) {
      
      if (tx_st != 1) {
        P2 <- if(tx_patient[i]) P_tx else P
      } else {
        P2 <- P
      }
      
      if (use_avg_start_dist & !is.null(avg_start_dist)) {
        
        diag_stage_i <- as.integer(site_data[i, "stage_imp"])
        month_i <- as.numeric(site_data[i, "month"])
        start_dist <- avg_start_dist %>% 
          filter(site==st & diagnosed_stage==diag_stage_i & month == month_i) %>% 
          select(X1, X2, X3, X4)
        
      } else {
        
        start_dist <- starting_distributions[[st]][i, -1]
      }
        
      # Run n_sim simulations for this patient
      for (sim in 1:n_sim) {
        diagnosed_stage <- simulate_single_patient(P2, start_dist, months=months) # remove id, only keep 4 stages
        
        if (!is.na(diagnosed_stage)) {
          # Accumulate for simulation-level results
          simulation_results[sim, diagnosed_stage] <- simulation_results[sim, diagnosed_stage] + 1
          
          # Accumulate for patient-level results
          patient_stage_counts[i, diagnosed_stage] <- patient_stage_counts[i, diagnosed_stage] + 1
          patient_diagnosed_counts[i] <- patient_diagnosed_counts[i] + 1
        }
      }
    }
    
    # Calculate simulation-level proportions (same as original code)
    for (sim in 1:n_sim) {
      n_diagnosed_in_sim <- sum(simulation_results[sim, ])
      if (n_diagnosed_in_sim > 0) {
        simulation_results[sim, ] <- simulation_results[sim, ] / n_diagnosed_in_sim
      } else {
        simulation_results[sim, ] <- rep(NA, 4)
      }
      simulation_diagnosis_rates[sim] <- n_diagnosed_in_sim / n_patients
    }
    
    # Calculate mean and confidence intervals across simulations (exactly like original)
    predicted_distributions[[st]] <- list(
      mean = colMeans(simulation_results, na.rm = TRUE),
      ci_lower = apply(simulation_results, 2, quantile, 0.025, na.rm = TRUE),
      ci_upper = apply(simulation_results, 2, quantile, 0.975, na.rm = TRUE),
      observed = as.numeric(table(factor(site_data$diagnosed_stage, levels = 1:4)) / nrow(site_data)),
      observed_lower = obs_ci(site_data)[["obs_lower"]],
      observed_upper = obs_ci(site_data)[["obs_upper"]],
      diagnosis_rate = mean(simulation_diagnosis_rates),
      diagnosis_rate_ci = quantile(simulation_diagnosis_rates, c(0.025, 0.975)),
      patient_level_results = list(
        e_patid = starting_distributions[[st]][, 1],
        stage_counts = patient_stage_counts,
        diagnosed_counts = patient_diagnosed_counts,
        stage_proportions = patient_stage_counts / pmax(patient_diagnosed_counts, 1)
      )
    )
  }
  
  return(predicted_distributions)
}


obs_ci <- function(site_data){
  
  library(boot)
  set.seed(123) 
  
  # Bootstrap function to calculate stage proportions
  stage_props <- function(data, indices) {
    # Resample the data
    resampled_stages <- data[indices]
    # Calculate proportions for each stage (1-4)
    props <- as.numeric(table(factor(resampled_stages, levels = 1:4)) / length(resampled_stages))
    return(props)
  }
  
  sites <- as.character(unique(site_data$site))
  
  boot_results <- boot(site_data$diagnosed_stage, stage_props, R = 10000)
  
  results_boot <- data.frame(
    stage = 1:4,
    percentage = as.numeric(table(factor(site_data$diagnosed_stage, levels = 1:4)) / nrow(site_data))
  )
  
  for(i in 1:4) {
    ci <- boot.ci(boot_results, type = "perc", index = i, conf = 0.95)
    results_boot$obs_lower[i] <- ci$percent[4]
    results_boot$obs_upper[i] <- ci$percent[5]
  }
  
  return(results_boot)
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
        observed_lower = site_data$observed_lower[stage],
        observed_upper = site_data$observed_upper[stage],
        difference = site_data$mean[stage] - site_data$observed[stage]
      ))
    }
  }
  
  # Create readable labels for cancer sites
  site_labels <- c(
    "panc" = "Pancreas",
    "oeso" = "Oesophagus",
    "stom" = "Stomach",
    "galb" = "Gallbladder"
  )
  comparison_df$cancer_site <- factor(comparison_df$cancer_site, 
                                      levels = names(site_labels),
                                      labels = site_labels)
  
  # Convert stage numbers to Roman numerals
  comparison_df$stage <- factor(comparison_df$stage,
                                levels = 1:4,
                                labels = c("I", "II", "III", "IV"))
  
  return(comparison_df)
}


# Visualization function
plot_comparison <- function(prediction_results, re_cali=FALSE) {
  
  library(ggplot2)
  
  if(!re_cali) comparison_df <- compare_distributions(prediction_results) else {
    
    comparison_df <- prediction_results
  }
  
  p <- ggplot(comparison_df, aes(x = factor(stage))) +
    geom_col(aes(y = observed, fill = "Observed"), alpha = 0.7, position = "dodge") +
    geom_col(aes(y = predicted, fill = "Predicted"), alpha = 0.7, position = "dodge2") +
    # Error bar for observed - positioned slightly to the left
    geom_errorbar(aes(ymin = observed_lower, ymax = observed_upper), 
                  width = 0.15, position = position_nudge(x = -0.08),
                  color = "#2E5F8A", linewidth = 0.6) +
    # Error bar for predicted - positioned slightly to the right
    geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper), 
                  width = 0.15, position = position_nudge(x = 0.08),
                  color = "#D9534F", linewidth = 0.6) +
    facet_wrap(~cancer_site, scales = "fixed") +
    labs(x = "Cancer Stage", y = "Proportion") +
    scale_fill_manual(name=NULL, values = c("Observed" = "steelblue", "Predicted" = "coral")) +
    ylim(0, 0.761) +
    theme_minimal()+
    theme(
      strip.text = element_text(size = 13, face = "bold"),  # Facet panel names
      axis.title.x = element_text(size = 13, face = "bold"),  # X axis label
      axis.title.y = element_text(size = 13, face = "bold"),  # Y axis label
      axis.text.x = element_text(size = 11),   # X axis tick labels
      axis.text.y = element_text(size = 11),   # Y axis tick labels
      legend.text = element_text(size = 12) 
    )
  
  return(p)
}

plot_prep_sub <- function(prediction_results, sub_name){
  
  sub_prediction_results <- list(
    stom = prediction_results[["stom"]][["sub_group"]][[sub_name]],
    panc = prediction_results[["panc"]][["sub_group"]][[sub_name]],
    oeso = prediction_results[["oeso"]][["sub_group"]][[sub_name]],
    galb = prediction_results[["galb"]][["sub_group"]][[sub_name]]
  )
  
  return(sub_prediction_results)
}


compare_distributions_interv2 <- function(pred_dist_interv) { 
  # for validation in high risk  
  comparison_df <- data.frame()
  
  for (site in names(pred_dist_interv)) {
    
    interv_data <- pred_dist_interv[[site]]
    
    for (stage in 1:4) {
      comparison_df <- rbind(comparison_df, data.frame(
        cancer_site = site,
        stage = stage,
        current = interv_data$observed[stage],
        current_lower = interv_data$observed_lower[stage],
        current_upper = interv_data$observed_upper[stage],
        intervention = interv_data$mean[stage],
        difference = interv_data$mean[stage] - interv_data$observed[stage],
        intervention_lower = interv_data$ci_lower[stage],
        intervention_upper = interv_data$ci_upper[stage]
      ))
    }
  }
  
  return(comparison_df)
}


plot_comparison_interv2 <- function(comparison_df) {
  # for validation in high risk
  
  library(ggplot2)
  
  # Create readable labels for cancer sites
  site_labels <- c(
    "panc" = "Pancreas",
    "oeso" = "Oesophagus",
    "stom" = "Stomach",
    "galb" = "Gallbladder"
  )
  comparison_df$cancer_site <- factor(comparison_df$cancer_site, 
                                      levels = names(site_labels),
                                      labels = site_labels)
  
  # Convert stage numbers to Roman numerals
  comparison_df$stage <- factor(comparison_df$stage,
                                levels = 1:4,
                                labels = c("I", "II", "III", "IV"))
  
  p <- ggplot(comparison_df, aes(x = factor(stage))) +
    geom_col(aes(y = current, fill = "Observed"), alpha = 0.9, position = "dodge") +
    geom_col(aes(y = intervention, fill = "Predicted"), alpha = 0.7, position = "dodge2") +
    geom_errorbar(aes(ymin = intervention_lower, ymax = intervention_upper),
                  width = 0.3, linewidth = 0.6, color = "#D9534F") +
    geom_errorbar(aes(ymin = current_lower, ymax = current_upper), 
                  width = 0.15, position = position_nudge(x = -0.08),
                  color = "#2E5F8A", linewidth = 0.6) +
    
    facet_wrap(~cancer_site, scales = "fixed") +
    labs(x = "Cancer Stage", y = "Proportion") +
    scale_fill_manual(name=NULL, values = c("Observed" = "steelblue", "Predicted" = "coral")) +
    ylim(0, 0.9) +
    theme_minimal()+
    theme(
      strip.text = element_text(size = 13, face = "bold"),  # Facet panel names
      axis.title.x = element_text(size = 13, face = "bold"),  # X axis label
      axis.title.y = element_text(size = 13, face = "bold"),  # Y axis label
      axis.text.x = element_text(size = 11),   # X axis tick labels
      axis.text.y = element_text(size = 11),   # Y axis tick labels
      legend.text = element_text(size = 13, face = "bold") 
    )
  
  return(p)
}



# Apply intervention effect -----------------------------------------------
compare_distributions_interv <- function(pred_dist_interv) {
  
  comparison_df <- data.frame()
  
  for (site in names(pred_dist_interv)) {

    interv_data <- pred_dist_interv[[site]]
    
    for (stage in 1:4) {
      comparison_df <- rbind(comparison_df, data.frame(
        cancer_site = site,
        stage = stage,
        current = interv_data$observed[stage],
        intervention = interv_data$mean[stage],
        difference = interv_data$mean[stage] - interv_data$observed[stage],
        intervention_lower = interv_data$ci_lower[stage],
        intervention_upper = interv_data$ci_upper[stage]
      ))
    }
  }
  
  return(comparison_df)
}


# Visualization function
plot_comparison_interv <- function(comparison_df) {
  library(ggplot2)
  
  # Create readable labels for cancer sites
  site_labels <- c(
    "panc" = "Pancreas",
    "oeso" = "Oesophagus",
    "stom" = "Stomach",
    "galb" = "Gallbladder"
  )
  comparison_df$cancer_site <- factor(comparison_df$cancer_site, 
                                      levels = names(site_labels),
                                      labels = site_labels)
  
  # Convert stage numbers to Roman numerals
  comparison_df$stage <- factor(comparison_df$stage,
                                levels = 1:4,
                                labels = c("I", "II", "III", "IV"))
  
  p <- ggplot(comparison_df, aes(x = factor(stage))) +
    geom_col(aes(y = current, fill = "Current"), alpha = 0.9, position = "dodge") +
    geom_col(aes(y = intervention, color = "Intervention"), 
             fill = NA, linewidth = 1, position = "dodge2") +
    geom_errorbar(aes(ymin = intervention_lower, ymax = intervention_upper),
                  width = 0.3, linewidth = 0.6, color = "#FF6B6B") +
    facet_wrap(~cancer_site, scales = "fixed") +
    labs(x = "Cancer Stage", y = "Proportion") +
    scale_fill_manual(name = "", 
                      values = c("Current" = "#008B8B"),
                      labels = c("Current")) +
    scale_color_manual(name = "",
                       values = c("Intervention" = "#FF6B6B"),
                       labels = c("Intervention")) +
    scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 2, override.aes = list(fill = NA, linewidth = 1))) +
    theme_minimal()+
    theme(
      strip.text = element_text(size = 13, face = "bold"),  # Facet panel names
      axis.title.x = element_text(size = 13, face = "bold"),  # X axis label
      axis.title.y = element_text(size = 13, face = "bold"),  # Y axis label
      axis.text.x = element_text(size = 11),   # X axis tick labels
      axis.text.y = element_text(size = 11),   # Y axis tick labels
      legend.text = element_text(size = 13, face = "bold") 
    )
  
  return(p)
}

# Validation old -----------------------------------------------------------

# This is not a good method to validate the model

# get the probability of correctly predicting the stage in the diagnosis month using optimized parameters
# or the probability of correctly predicting the undiagnosed stage 1 month before diagnosis, i.e. undiagnosed at same cancer stage
predict_stage_at_diagnosis <- function(optimized_params, month_diagnosis, 
                                       diagnosed_stage, before_diag = FALSE, 
                                       n_cycles = 24, n_sim = 10000) {
  # Extract parameters
  tp12 <- optimized_params[1]
  tp23 <- optimized_params[2]
  tp34 <- optimized_params[3]
  tp1 <- optimized_params[4]
  tp2 <- optimized_params[5]
  tp3 <- optimized_params[6]
  tp4 <- optimized_params[7]
  
  # Create transition matrix
  P <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
  
  # Simulate for each possible starting stage
  stage_predictions <- rep(0, 4)  # 4 possible starting_stages
  
  for (start_stage in 1:4) {
    # Simulate Markov chain
    sim_results <- simulate_markov(P, initial_state = start_stage, n_cycles = n_cycles)
    
    # At the diagnosis month, get probability of being diagnosed at each stage
    if (!before_diag) {
      # the month of diagnosis, the state is diagnosed
      stage_predictions[start_stage] <- sim_results[month_diagnosis + 1, diagnosed_stage + 4] # diagnosed states
    } else {
      # 1 month before the diagnosis, so the state is undiagnosed and 
      # cancer stage assumed to be the same
      stage_predictions[start_stage] <- sim_results[month_diagnosis, diagnosed_stage] # undiagnosed states
    }
  }
  
  # Calculate weighted overall stage distribution for each diagnosed stage
  # Calculate the proper initial stage distribution for this diagnosed stage
  stage_weights <- initial_stage_distribution(month_diagnosis, diagnosed_stage, 
                                                tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                                n_sim = n_sim)
  
  # Weight the predictions by the initial stage distribution
  overall_predictions<- sum(stage_predictions* stage_weights)
  
  return(overall_predictions)
}


# Modified validation function to find the state with largest probability
# Most simplified version - returns only the state index (1-8)
predict_most_likely_state_index <- function(optimized_params, month_diagnosis, 
                                            diagnosed_stage, before_diag = FALSE, 
                                            n_cycles = 24, n_sim = 10000) {
  # Extract parameters
  tp12 <- optimized_params[1]
  tp23 <- optimized_params[2]
  tp34 <- optimized_params[3]
  tp1 <- optimized_params[4]
  tp2 <- optimized_params[5]
  tp3 <- optimized_params[6]
  tp4 <- optimized_params[7]
  
  # Create transition matrix
  P <- create_transition_matrix(tp12, tp23, tp34, tp1, tp2, tp3, tp4)
  
  # Calculate the proper initial stage distribution for this diagnosed stage
  stage_weights <- initial_stage_distribution(month_diagnosis, diagnosed_stage, 
                                                tp12, tp23, tp34, tp1, tp2, tp3, tp4,
                                                n_sim = n_sim)
  
  # Initialize overall state probabilities (8 states total)
  overall_state_probs <- rep(0, 8)
  
  # Simulate for each possible starting stage and weight by initial distribution
  for (start_stage in 1:4) {
    # Only process if this starting stage has non-zero weight
    if (stage_weights[start_stage] > 0) {
      # Simulate Markov chain
      sim_results <- simulate_markov(P, initial_state = start_stage, n_cycles = n_cycles)
      
      # Get state probabilities at the target month
      if (!before_diag) {
        # At the diagnosis month
        target_month_probs <- sim_results[month_diagnosis + 1, ]
      } else {
        # 1 month before diagnosis
        target_month_probs <- sim_results[month_diagnosis, ]
      }
      
      # Add weighted probabilities to overall distribution
      overall_state_probs <- overall_state_probs + (stage_weights[start_stage] * target_month_probs)
    }
  }
  
  # Return the state index with maximum probability
  return(which.max(overall_state_probs))
}


# Function to create comprehensive validation dot plot with separate panels
create_validation_dotplot <- function(d4s, title_suffix = "") {
  
  # Prepare data for plotting - sort patients by diagnosed stage within each site
  plot_data <- d4s %>%
    arrange(site, diagnosed_stage, e_patid) %>%
    group_by(site) %>%
    mutate(
      # Create patient order within each site
      patient_order = row_number(),
      # Convert diagnosed_stage to reference state (diagnosed stages 5-8)
      reference_state = diagnosed_stage + 4,
      # Reference for 1 month before should be undiagnosed version
      reference_1monthb4 = diagnosed_stage
    ) %>%
    ungroup()
  
  # Create the main dot plot with facets
  p1 <- ggplot(plot_data, aes(x = patient_order)) +
    # Reference line (actual diagnosed state)
    geom_line(aes(y = reference_state), color = "red", linewidth = 1, alpha = 0.8) +
    # Reference line for 1 month before (undiagnosed state)
    geom_line(aes(y = reference_1monthb4), color = "blue", linewidth = 1, alpha = 0.8) +
    # Error bars for predictions
    geom_errorbar(aes(ymin = most_likely_state, ymax = most_likely_state), 
                  color = "darkgreen", width = 0.5, alpha = 0.7) +
    geom_errorbar(aes(ymin = most_likely_state_1monthb4, ymax = most_likely_state_1monthb4), 
                  color = "orange", width = 0.5, alpha = 0.7) +
    # Points for predictions
    geom_point(aes(y = most_likely_state), color = "darkgreen", size = 1.5, alpha = 0.8) +
    geom_point(aes(y = most_likely_state_1monthb4), color = "orange", size = 1.5, alpha = 0.8) +
    # Separate panels for each cancer site
    facet_wrap(~ site, ncol = 1, scales = "free_x", 
               labeller = labeller(site = function(x) paste("Cancer Site:", toupper(x)))) +
    scale_y_continuous(
      breaks = 1:8,
      labels = c("Undiag S1", "Undiag S2", "Undiag S3", "Undiag S4",
                 "Diag S1", "Diag S2", "Diag S3", "Diag S4"),
      limits = c(0.5, 8.5)
    ) +
    labs(
      title = paste("Validation: Predicted vs Reference States by Cancer Site", title_suffix),
      subtitle = "Patients sorted by diagnosed stage within each site",
      x = "Patient Order (sorted by diagnosed stage)",
      y = "State (1-4: Undiagnosed, 5-8: Diagnosed)",
      caption = "Red line: Reference (actual diagnosed stage). Blue line: Reference 1 month before (undiagnosed).\nGreen points: Prediction at diagnosis. Orange points: Prediction 1 month before."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  
  return(p1)
}

# Function to create accuracy summary plot
create_accuracy_summary <- function(d4s, title_suffix = "") {
  
  # Calculate accuracy metrics
  accuracy_data <- d4s %>%
    mutate(
      reference_state = diagnosed_stage + 4,
      # Check if predictions match reference
      at_diagnosis_exact = (most_likely_state == reference_state),
      one_month_before_exact = (most_likely_state_1monthb4 == diagnosed_stage), # compare to undiagnosed stage
      # Check if predictions get the stage right (regardless of diagnosis status)
      at_diagnosis_stage = ((most_likely_state - 1) %% 4 + 1) == diagnosed_stage,
      one_month_before_stage = ((most_likely_state_1monthb4 - 1) %% 4 + 1) == diagnosed_stage
    ) %>%
    group_by(site) %>%
    summarise(
      n_patients = n(),
      exact_match_at_diag = mean(at_diagnosis_exact) * 100,
      exact_match_1month_before = mean(one_month_before_exact) * 100,
      stage_correct_at_diag = mean(at_diagnosis_stage) * 100,
      stage_correct_1month_before = mean(one_month_before_stage) * 100,
      .groups = 'drop'
    ) %>%
    pivot_longer(
      cols = contains("_correct") | contains("match"),
      names_to = "metric",
      values_to = "accuracy_pct"
    ) %>%
    mutate(
      metric = factor(metric, 
                      levels = c("exact_match_at_diag", "exact_match_1month_before", 
                                 "stage_correct_at_diag", "stage_correct_1month_before"),
                      labels = c("Exact Match\n(At Diagnosis)", "Exact Match\n(1 Month Before)",
                                 "Stage Correct\n(At Diagnosis)", "Stage Correct\n(1 Month Before)"))
    )
  
  # Create accuracy plot
  p2 <- ggplot(accuracy_data, aes(x = metric, y = accuracy_pct, fill = site)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(round(accuracy_pct, 1), "%")), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3) +
    scale_fill_brewer(type = "qual", palette = "Set2", name = "Cancer Site") +
    scale_y_continuous(labels = percent_format(scale = 1), limits = c(0, 100)) +
    labs(
      title = paste("Prediction Accuracy by Cancer Site", title_suffix),
      x = "Accuracy Metric",
      y = "Accuracy (%)",
      caption = "Exact Match: Predicted state exactly matches reference (including diagnosis status).\nStage Correct: Predicted cancer stage matches actual stage (regardless of diagnosis status)."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  return(p2)
}

# Function to create stage-specific comparison
create_stage_comparison <- function(d4s, title_suffix = "") {
  
  # Prepare data for stage-specific analysis
  stage_data <- d4s %>%
    mutate(
      reference_state = diagnosed_stage + 4,
      # Create stage difference metrics
      diff_at_diag = most_likely_state - reference_state,
      diff_1month_before = most_likely_state_1monthb4 - diagnosed_stage
    ) %>%
    pivot_longer(
      cols = c(diff_at_diag, diff_1month_before),
      names_to = "prediction_type",
      values_to = "difference"
    ) %>%
    mutate(
      prediction_type = factor(prediction_type,
                               levels = c("diff_at_diag", "diff_1month_before"),
                               labels = c("At Diagnosis", "1 Month Before")),
      difference_category = case_when(
        difference == 0 ~ "Perfect Match",
        difference > 0 ~ "Over-predicted",
        difference < 0 ~ "Under-predicted"
      )
    )
  
  # Create stage comparison plot
  p3 <- ggplot(stage_data, aes(x = factor(diagnosed_stage), fill = difference_category)) +
    geom_bar(position = "fill", alpha = 0.8) +
    facet_grid(prediction_type ~ site, scales = "free_x") +
    scale_fill_manual(
      values = c("Perfect Match" = "darkgreen", "Over-predicted" = "orange", "Under-predicted" = "red"),
      name = "Prediction Quality"
    ) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title = paste("Prediction Quality by Diagnosed Stage and Site", title_suffix),
      x = "Diagnosed Stage",
      y = "Proportion of Patients",
      caption = "Perfect Match: Predicted exactly right.\nOver-predicted: Predicted higher stage than actual.\nUnder-predicted: Predicted lower stage than actual."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  return(p3)
}


