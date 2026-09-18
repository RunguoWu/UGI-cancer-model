# PSA -----
# Intialisation -----

rm(list = ls())

source("/data/WIPH-CanDetect/HealthEco/Martin/route.R")
source(file.path(scr, "tx_model", "fn_analysis.R"))
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))


# Load libraries

library(future)
library(future.apply)
library(tidyverse)
library(survival)
library(gtsummary)
library(skimr)
library(dplyr)
library(rstpm2)
library(officer)
library(parallelly)

# New functions

predict_stage_distribution3 <- function(patient_data, params, n_sim = 500, months = 24,
                                        tx, tx_prob = 1, use_avg_start_dist = FALSE, 
                                        avg_start_dist = NULL) {
  
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
      scaled_tp4 <- 1 - (1 - optimized_params[4])^tx_st
      scaled_tp5 <- 1 - (1 - optimized_params[5])^tx_st
      scaled_tp6 <- 1 - (1 - optimized_params[6])^tx_st
      
      P_tx <- create_transition_matrix(
        optimized_params[1], optimized_params[2], optimized_params[3],
        scaled_tp4,          scaled_tp5,          scaled_tp6, 
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

extract_tx <- function(model) {
  tx_img <- 1/exp(coef(model)["ng12_red_flag2Imaging"])
  tx_2ww <- 1/exp(coef(model)["ng12_red_flag22 Week Wait"])
  rt <- c(tx_img, tx_2ww)
  names(rt) <- c("imaging", "2ww")
  return(rt)
}


# Base case ----

# 1. Load data, set names and formulas

pop_nam   <- "symp_stageImputed_upd2026"
char_name <- c("age70plus", "female")
tx_name   <- c("2ww", "imaging")
site_name <- c("panc", "oeso", "stom", "galb")

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds")))

d4s_no_red <- d4s %>%
  filter(ng12_red_flag == "No red flag")

d4s_surv <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) %>% 
  mutate(death = if_else(is.na(death_date), 0, 1),
         death_cancer = if_else(is.na(death_cancer), 0, death_cancer),
         death_upGI = if_else(is.na(death_upGI), 0, death_cancer),
         death_oeso = if_else(is.na(death_oeso), 0, death_cancer),
         death_stom = if_else(is.na(death_stom), 0, death_cancer),
         death_panc = if_else(is.na(death_panc), 0, death_cancer),
         death_galb = if_else(is.na(death_galb), 0, death_cancer)) %>% 
  mutate(fu_diag2cens = as.numeric(fu_end_date - cancerdate),
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 0.5, fu_diag2cens),
         fu_diag2cens = fu_diag2cens/365.25) %>% 
  select(e_patid, site, fu_diag2cens, death_cancer, age10_cent60, nonwhite, female, imd5_imp2, stage_imp)

fm_tx <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")

# 2. Survival model

sx_mod_list <- list()

for (st in site_name) {
  
  d4s_sub <- subset(d4s_surv, site == st)
  
  if (st == "panc") {
    fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 * stage_imp + nonwhite + female + imd5_imp2 + stage_imp
  } else {
    fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 + nonwhite + female + imd5_imp2 + stage_imp
  }
  
  fit               <- stpm2(fm, d4s_sub, df = 5)
  sx_mod_list[[st]] <- fit
  
}

# 3. Treatment-effect model (by site only)

tx_mod_list <- list()
tx_list <- list()

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    for (st in site_name) {
      
      d4s_sub        <- subset(d4s, site == st)
      sub_name       <- paste0("age70plus", char_value1, "_", "female", char_value2)
      time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
      wei            <- survreg(fm_tx, d4s_sub, dist = "weibull")
      
      tx_mod_list[[sub_name]][[st]] <- wei
      tx_list[[sub_name]][[st]]     <- extract_tx(wei)
      
    }
  }
}


# 4. Transition probabilities

params_list <- readRDS(file.path(wd, "params_list_upd2026.rds"))


# 5. Predict stage shift 

tx_name <- c("2ww", "imaging")

output_list <- list()
ind_output_list <- list()

for (tx_strat in tx_name) {
  ind_output_list[[tx_strat]] <- list()
  
  for (char_value1 in 0:1) {
    for (char_value2 in 0:1){
      
      sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
      
      sub_data <- d4s_no_red %>%
        filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
      
      params <- params_list[[sub_name]]
      
      tx <- sapply(tx_list[[sub_name]], "[", tx_strat)
      names(tx) <- gsub(paste0("\\.", tx_strat, "$"), "", names(tx))
      
      set.seed(123)
      pred_dist_interv <- predict_stage_distribution3(sub_data, params,  n_sim = 1000, months = 24, tx = tx, tx_prob = 1)
      comparison_df    <- compare_distributions_interv(pred_dist_interv)
      
      output_list[[tx_strat]][[sub_name]]     <- comparison_df
      ind_output_list[[tx_strat]][[sub_name]] <- pred_dist_interv
    }
  }
}

ind_output_list_img <- ind_output_list[["imaging"]]
ind_output_list_2ww <- ind_output_list[["2ww"]]

stage_list <- list()

for (st in site_name){
  
  stage_list[[st]] <- NULL
  
  for (i in 1:length(ind_output_list_img)) {
    
    e_patid <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["e_patid"]]
    
    stage_proportions_img <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
    stage_proportions_2ww <- ind_output_list_2ww[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
    
    stage_proportions_df <- data.frame(
      e_patid = e_patid,
      X1_img = stage_proportions_img[, 1],
      X2_img = stage_proportions_img[, 2],
      X3_img = stage_proportions_img[, 3],
      X4_img = stage_proportions_img[, 4],
      X1_2ww = stage_proportions_2ww[, 1],
      X2_2ww = stage_proportions_2ww[, 2],
      X3_2ww = stage_proportions_2ww[, 3],
      X4_2ww = stage_proportions_2ww[, 4],
      stringsAsFactors = FALSE)
    
    stage_list[[st]] <- rbind(stage_list[[st]], stage_proportions_df)
    
  }
}

stage_df <- do.call(rbind, stage_list)

d4s2 <- d4s_no_red %>%
  select(e_patid, site, female, age10_cent60, age10_new, age70plus, nonwhite, imd5_imp2, stage_imp, death_cancer) %>%
  left_join(stage_df)

# 5. Predict survival

surv_years <- 5

pred_list <- list()

for (st in site_name) {
  
  mod     <- sx_mod_list[[st]]
  st_data <- d4s2 %>% filter(site == st)
  
  pred_list[[st]] <- list()
  pred_st <- NULL
  
  for (stage in as.character(1:4)) {
    
    newdata_int <- st_data %>%
      mutate(fu_diag2cens = surv_years, stage_imp = stage)
    
    pred_int <- predict(mod, newdata = newdata_int, type = "surv")
    pred_st  <- cbind(pred_st, pred_int)
    
  }
  
  newdata_cur <- st_data %>% mutate(fu_diag2cens = surv_years)
  pred_cur    <- predict(mod, newdata = newdata_cur, type = "surv")
  
  surv_st           <- data.frame(st_data$e_patid, pred_st, pred_cur)
  colnames(surv_st) <- c("e_patid", paste0("diag_s", 1:4), "surv_cur")
  
  pred_list[[st]] <- surv_st
  
}

pred_dt <- do.call(rbind, pred_list)

d4s2 <- d4s2 %>%
  left_join(pred_dt)

# Run scenarios

s_img_all     <- partial_intervention_summary(d4s2, pct_img = 1, pct_2ww = 0, seed = 123)
s_2ww_all     <- partial_intervention_summary(d4s2, pct_img = 0, pct_2ww = 1, seed = 123)
s_img_haf     <- partial_intervention_summary(d4s2, pct_img = 0.5, pct_2ww = 0, seed = 123)
s_2ww_haf     <- partial_intervention_summary(d4s2, pct_img = 0, pct_2ww = 0.5, seed = 123)
s_img_2ww_qua <- partial_intervention_summary(d4s2, pct_img = 0.25, pct_2ww = 0.25, seed = 123)

rt <- cbind(s_img_all[, c("site", "stage", "n_total", "obs_stage", "obs_surv", "pred_stage")],
            s_img_all[, c("pred_surv")],
            s_2ww_all[, c("pred_stage")],
            s_2ww_all[, c("pred_surv")],
            s_img_haf[, c("pred_stage")],
            s_img_haf[, c("pred_surv")],
            s_2ww_haf[, c("pred_stage")],
            s_2ww_haf[, c("pred_surv")],
            s_img_2ww_qua[, c("pred_stage")],
            s_img_2ww_qua[, c("pred_surv")])

rt[, 4:15] <- round(rt[, 4:15], 3) * 100

colnames(rt) <- c("Site", "Stage", "N", "Stage distribution", "Five-year survival",
                  "pred_stage_img_all",
                  "pred_surv_img_all",
                  "pred_stage_2ww_all",
                  "pred_surv_2ww_all",
                  "pred_stage_img_haf",
                  "pred_surv_img_haf",
                  "pred_stage_2ww_haf",
                  "pred_surv_2ww_haf",
                  "pred_stage_img_2ww_qua",
                  "pred_surv_img_2ww_qua")

print(rt)

write.csv(rt, file.path(output, "MV_stagedist_5yrsurv_2026-07.csv"))

rm(mod_list, ind_output_list_img, ind_output_list_2ww, tx_list,
   stage_list, stage_df, pred_list, pred_dt, d4s2, s_img_all,
   s_2ww_all, s_img_haf, s_2ww_haf, s_img_2ww_qua, rt)


# Set up the PSA ----
# Set the configurations

pop_nam   <- "symp_stageImputed_upd2026"
char_name <- c("age70plus", "female")
tx_name   <- c("2ww", "imaging")
site_name <- c("panc", "oeso", "stom", "galb")

surv_years <-5

df_choices <- list(oeso = 5, stom = 5, panc = 5, galb = 5)

fm_tx <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")


## Transition probabilities
#
# dt_list <- list()
# for (site in site_name) {
# 
#   dt_list[[paste0(site, "_", char_name[1], 0, "_", char_name[2], 0)]] <-
#     readRDS(file.path("/data/WIPH-CanDetect/HealthEco/output", "optim",
#                       paste0(opt_tag, "_",
#                              char_name[1], 0, "_",
#                              char_name[2], 0, "_", site, ".rds")))
# 
#   dt_list[[paste0(site, "_", char_name[1], 0, "_", char_name[2], 1)]] <-
#     readRDS(file.path("/data/WIPH-CanDetect/HealthEco/output", "optim",
#                       paste0(opt_tag, "_",
#                              char_name[1], 0, "_",
#                              char_name[2], 1, "_", site, ".rds")))
# 
#   dt_list[[paste0(site, "_", char_name[1], 1, "_", char_name[2], 0)]] <-
#     readRDS(file.path("/data/WIPH-CanDetect/HealthEco/output", "optim",
#                       paste0(opt_tag, "_",
#                              char_name[1], 1, "_",
#                              char_name[2], 0, "_", site, ".rds")))
# 
#   dt_list[[paste0(site, "_", char_name[1], 1, "_", char_name[2], 1)]] <-
#     readRDS(file.path("/data/WIPH-CanDetect/HealthEco/output", "optim",
#                       paste0(opt_tag, "_",
#                              char_name[1], 1, "_",
#                              char_name[2], 1, "_", site, ".rds")))
# }
# 
# total_runs <- length(dt_list[[1]]$params_record)
# global_keep_idx <- rep(TRUE, total_runs)
# 
# for (nm in names(dt_list)) {
#   all_records <- dt_list[[nm]]$params_record
#   
#   for (i in 1:total_runs) {
#     
#     if (!global_keep_idx[i]) next
#     
#     p <- all_records[[i]]$par_cur
#     tp12 <- p[1]; tp23 <- p[2]; tp34 <- p[3]
#     tp1  <- p[4]; tp2  <- p[5]; tp3  <- p[6]; tp4  <- p[7]
#     
#     v_comp  <- (tp12 + tp1 <= 1) && (tp23 + tp2 <= 1) && (tp34 + tp3 <= 1) && (tp4 <= 1)
#     v_accum <- (tp1 <= tp2) && (tp2 <= tp3) && (tp3 <= tp4)
#     v_progr <- (tp12 <= tp23) && (tp23 <= tp34)
#     v_range <- all(p >= 0.01) && all(p <= 1)
#     
#     if (!(v_comp && v_accum && v_progr && v_range)) {
#       global_keep_idx[i] <- FALSE
#     }
#   }
# }
# 
# for (nm in names(dt_list)) {
#   dt_list[[nm]]$params_record <- dt_list[[nm]]$params_record[global_keep_idx]
# }
# 
# valid_indices_lookup <- list()
# 
# for (nm in names(dt_list)) {
#   llh_records <- sapply(dt_list[[nm]]$params_record, "[[", "value_cur")
#   threshold <- quantile(llh_records, probs = 0.01)
#   valid_indices_lookup[[nm]] <- which(llh_records <= threshold)
# }
# 
# n_psa_tp_samples <- 1000
# params_all <- vector("list", n_psa_tp_samples)
# 
# set.seed(20260701)
# 
# for (idx in 1:n_psa_tp_samples) {
#   params_list_idx <- list()
#   
#   for (char_value1 in 0:1) {
#     for (char_value2 in 0:1) {
#       sub_name <- paste0(char_name[1], char_value1, "_", char_name[2], char_value2)
#       keys <- names(dt_list)[grepl(sub_name, names(dt_list))]
#       
#       params <- do.call(rbind, lapply(keys, function(nm) {
#         
#         valid_indices <- valid_indices_lookup[[nm]]
#         
#         random_safe_idx <- sample(valid_indices, 1)
#         selected_record <- dt_list[[nm]][["params_record"]][[random_safe_idx]]
#         
#         data.frame(site = sub("_.*", "", nm), t(selected_record$par_cur))
#       }))
#       
#       colnames(params) <- c("site", "tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
#       params_list_idx[[sub_name]] <- params
#     }
#   }
#   params_all[[idx]] <- params_list_idx
# }
# 
# rm(dt_list)
# 
# saveRDS(params_all, file.path(output, paste0("params_all_for_bootstrap.rds")))

params_all <- readRDS(file.path(output, paste0("params_all_for_bootstrap.rds")))



# PSA function ---- 

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds")))

psa <- function(b, d4s, params_all, df_choices, fm_tx, site_name, surv_years) {
  
  tp_log_b <- list()   
  
  result <- tryCatch({
    
    #### 1. Bootstrap ####
    
    d4s <- d4s %>%
      group_by(site) %>%
      sample_frac(size = 1, replace = TRUE) %>%
      ungroup() %>%
      mutate(e_patid = paste0(e_patid, "_", row_number()))
    
    boot_d4s <- d4s 
    
    boot_d4s_no_red <- boot_d4s %>% 
      filter(ng12_red_flag == "No red flag")
    
    boot_d4s_surv <- boot_d4s %>%
      mutate(death = if_else(is.na(death_date), 0, 1),
             death_cancer = if_else(is.na(death_cancer), 0, death_cancer),
             death_upGI = if_else(is.na(death_upGI), 0, death_cancer),
             death_oeso = if_else(is.na(death_oeso), 0, death_cancer),
             death_stom = if_else(is.na(death_stom), 0, death_cancer),
             death_panc = if_else(is.na(death_panc), 0, death_cancer),
             death_galb = if_else(is.na(death_galb), 0, death_cancer)) %>% 
      mutate(fu_diag2cens = as.numeric(fu_end_date - cancerdate),
             fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 0.5, fu_diag2cens),
             fu_diag2cens = fu_diag2cens/365.25) %>% 
      select(e_patid, site, fu_diag2cens, death_cancer, age10_cent60, nonwhite, female, imd5_imp2, stage_imp)
    
    #### 2. Survival model ####
    
    sx_mod_list <- list()
    
    for (st in site_name) {
      
      d4s_sub <- subset(boot_d4s_surv, site == st)
      
      if (st == "panc") {
        fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 * stage_imp + nonwhite + female + imd5_imp2 + stage_imp
      } else {
        fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 + nonwhite + female + imd5_imp2 + stage_imp
      }
      
      fit               <- stpm2(fm, d4s_sub, df = 5)
      sx_mod_list[[st]] <- fit
      
    }
    
    #### 3. Treatment-effect model (by site only) #### 
    
    tx_mod_list <- list()
    tx_list <- list()
    
    for (char_value1 in 0:1) {
      for (char_value2 in 0:1){
        for (st in site_name) {
          
          d4s_sub        <- subset(boot_d4s, site == st)
          sub_name       <- paste0("age70plus", char_value1, "_", "female", char_value2)
          d4s_sub$time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
          wei            <- survreg(fm_tx, d4s_sub, dist = "weibull")
          
          tx_mod_list[[sub_name]][[st]] <- wei
          tx_list[[sub_name]][[st]]     <- extract_tx(wei)
          
        }
      }
    }
    
    #### 3. Transition probabilities#### 
    
    selected_index <- sample(1:1000, size = 1)
    params_list <- params_all[[selected_index]]
    
    # params_list <- params_list
    
    #### 4. Predict stage shift ####
    
    output_list <- list()
    ind_output_list <- list()
    
    for (tx_strat in tx_name) {
      ind_output_list[[tx_strat]] <- list()
      
      for (char_value1 in 0:1) {
        for (char_value2 in 0:1){
          
          sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
          
          sub_data <- boot_d4s_no_red %>%
            filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
          
          params <- params_list[[sub_name]]
          
          tx <- sapply(tx_list[[sub_name]], "[", tx_strat)
          names(tx) <- gsub(paste0("\\.", tx_strat, "$"), "", names(tx))
          
          pred_dist_interv <- predict_stage_distribution3(sub_data, params,  n_sim = 1000, months = 24, tx = tx, tx_prob = 1)
          comparison_df    <- compare_distributions_interv(pred_dist_interv)
          
          output_list[[tx_strat]][[sub_name]]     <- comparison_df
          ind_output_list[[tx_strat]][[sub_name]] <- pred_dist_interv
        }
      }
    }
    
    ind_output_list_img <- ind_output_list[["imaging"]]
    ind_output_list_2ww <- ind_output_list[["2ww"]]
    
    stage_list <- list()
    
    for (st in site_name){
      
      stage_list[[st]] <- NULL
      
      for (i in 1:length(ind_output_list_img)) {
        
        e_patid <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["e_patid"]]
        
        stage_proportions_img <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
        stage_proportions_2ww <- ind_output_list_2ww[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
        
        stage_proportions_df <- data.frame(
          e_patid = e_patid,
          X1_img = stage_proportions_img[, 1],
          X2_img = stage_proportions_img[, 2],
          X3_img = stage_proportions_img[, 3],
          X4_img = stage_proportions_img[, 4],
          X1_2ww = stage_proportions_2ww[, 1],
          X2_2ww = stage_proportions_2ww[, 2],
          X3_2ww = stage_proportions_2ww[, 3],
          X4_2ww = stage_proportions_2ww[, 4],
          stringsAsFactors = FALSE)
        
        stage_list[[st]] <- rbind(stage_list[[st]], stage_proportions_df)
        
      }
    }
    
    stage_df <- do.call(rbind, stage_list)
    
    boot_d4s2 <- boot_d4s_no_red %>%
      select(e_patid, site, female, age10_cent60, age10_new, age70plus, nonwhite, imd5_imp2, stage_imp, death_cancer) %>%
      left_join(stage_df)
    
    
    #### 5. Predict survival ####
    
    pred_list <- list()
    
    for (st in site_name) {
      
      mod     <- sx_mod_list[[st]]
      st_data <- boot_d4s2 %>% filter(site == st)
      
      pred_list[[st]] <- list()
      pred_st <- NULL
      
      for (stage in as.character(1:4)) {
        
        newdata_int <- st_data %>%
          mutate(fu_diag2cens = surv_years, stage_imp = stage)
        
        pred_int <- predict(mod, newdata = newdata_int, type = "surv")
        pred_st  <- cbind(pred_st, pred_int)
        
      }
      
      newdata_cur <- st_data %>% mutate(fu_diag2cens = surv_years)
      pred_cur    <- predict(mod, newdata = newdata_cur, type = "surv")
      
      surv_st           <- data.frame(st_data$e_patid, pred_st, pred_cur)
      colnames(surv_st) <- c("e_patid", paste0("diag_s", 1:4), "surv_cur")
      
      pred_list[[st]] <- surv_st
      
    }
    
    pred_dt <- do.call(rbind, pred_list)
    
    boot_d4s2 <- boot_d4s2 %>%
      left_join(pred_dt)
    
    #### 6. Run scenarios ####
    
    s_img_all     <- partial_intervention_summary(boot_d4s2, pct_img = 1, pct_2ww = 0, seed = b)
    s_2ww_all     <- partial_intervention_summary(boot_d4s2, pct_img = 0, pct_2ww = 1, seed = b)
    s_img_haf     <- partial_intervention_summary(boot_d4s2, pct_img = 0.5, pct_2ww = 0, seed = b)
    s_2ww_haf     <- partial_intervention_summary(boot_d4s2, pct_img = 0, pct_2ww = 0.5, seed = b)
    s_img_2ww_qua <- partial_intervention_summary(boot_d4s2, pct_img = 0.25, pct_2ww = 0.25, seed = b)
    
    #### 7. Clean up ####
    
    rt_b <- cbind(s_img_all[, c("site", "stage", "n_total", "obs_stage", "obs_surv", "pred_stage")], 
                  s_img_all[, c("pred_surv")],
                  s_2ww_all[, c("pred_stage")],
                  s_2ww_all[, c("pred_surv")],
                  s_img_haf[, c("pred_stage")],
                  s_img_haf[, c("pred_surv")],
                  s_2ww_haf[, c("pred_stage")],
                  s_2ww_haf[, c("pred_surv")],
                  s_img_2ww_qua[, c("pred_stage")],
                  s_img_2ww_qua[, c("pred_surv")]
    )
    
    rt_b[, 4:15] <- round(rt_b[, 4:15], 3) * 100
    colnames(rt_b) <- c("Site", "Stage", "N", "Stage distribution", "Five-year survival", 
                        "pred_stage_img_all",  
                        "pred_surv_img_all", 
                        "pred_stage_2ww_all",
                        "pred_surv_2ww_all", 
                        "pred_stage_img_haf", 
                        "pred_surv_img_haf",
                        "pred_stage_2ww_haf",
                        "pred_surv_2ww_haf",
                        "pred_stage_img_2ww_qua", 
                        "pred_surv_img_2ww_qua"
    )
    
    rt_b <- cbind(psa_run = b, rt_b)
    
    list(status = "ok", rt = rt_b, error_msg = NA_character_)
    
  }, error = function(e) {
    message(sprintf("psa_run %d failed and was skipped: %s", b, conditionMessage(e)))
    list(status = "failed", rt = NULL, error_msg = conditionMessage(e))
  })
  
  # tp_log_b is populated regardless of success/failure, since it lives in
  # psa()'s own frame, untouched by which branch of tryCatch executed.
  tp_log_out <- if (length(tp_log_b) > 0) do.call(rbind, tp_log_b) else NULL
  
  c(result, list(tp_log = tp_log_out, psa_run = b))
}

# Run PSA -----

num_cores <- parallelly::availableCores() -1
message(sprintf("Using %d worker(s) for future_lapply (plan: multisession).", num_cores))

plan(multisession, workers = num_cores)

# Run the PSA

t1 <- Sys.time()
t1

n_psa <- 1000  

psa_results_list <- future_lapply(
  1:n_psa,
  function(b) {
    psa(b = b, d4s = d4s, params_all = params_all, df_choices = df_choices, 
        fm_tx = fm_tx, site_name = site_name, surv_years = surv_years
    )
  }, future.seed = TRUE
)

t2 <- Sys.time()
print(t2 - t1)

# Reset parallel processing

plan(sequential)

# Drop skipped (NULL) iterations and report how many failed

n_failed <- sum(sapply(psa_results_list, is.null))
cat(sprintf("%d of %d PSA iterations failed and were skipped.\n", n_failed, n_psa))

psa_results_list <- Filter(Negate(is.null), psa_results_list)

psa_master  <- do.call(rbind, lapply(psa_results_list, `[[`, "rt"))
tp_log_master <- do.call(rbind, lapply(psa_results_list, `[[`, "tp_log"))

saveRDS(psa_master, file.path(output, "psa_results_upd2026.rds"))


# Inspect the tp log ----

# Overall breach rate by site/scenario

status_summary <- data.frame(
  psa_run = sapply(psa_results_list, `[[`, "psa_run"),
  status  = sapply(psa_results_list, `[[`, "status")
)

status_summary


# Summarise ----

stage_cols <- c("Stage distribution", grep("^pred_stage_", colnames(psa_master), value = TRUE))

stage_summary <- psa_master %>%
  group_by(Site, Stage) %>%
  summarise(across(all_of(stage_cols),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        lower = ~ quantile(.x, probs = 0.025, na.rm = TRUE),
                        upper = ~ quantile(.x, probs = 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),.groups = "drop")

panc_stage <- stage_summary %>%
  filter(Site == "panc") %>%
  transmute(
    Site, Stage,
    `Stage distribution`      = sprintf("%.1f (%.1f-%.1f)", `Stage distribution_mean`, `Stage distribution_lower`, `Stage distribution_upper`),
    `Image (All)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_all_mean, pred_stage_img_all_lower, pred_stage_img_all_upper),
    `2WW (All)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_all_mean, pred_stage_2ww_all_lower, pred_stage_2ww_all_upper),
    `Image (50%)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_haf_mean, pred_stage_img_haf_lower, pred_stage_img_haf_upper),
    `2WW (50%)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_haf_mean, pred_stage_2ww_haf_lower, pred_stage_2ww_haf_upper),
    `Image (25%) + 2WW (25%)`  = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_2ww_qua_mean, pred_stage_img_2ww_qua_lower, pred_stage_img_2ww_qua_upper)
  ) %>%
  print(width = Inf)

write.csv(panc_stage, file.path(output, "MV_stagedist_panc_stage_2026-09.csv"))


oeso_stage <- stage_summary %>%
  filter(Site == "oeso") %>%
  transmute(
    Site, Stage,
    `Stage distribution`      = sprintf("%.1f (%.1f-%.1f)", `Stage distribution_mean`, `Stage distribution_lower`, `Stage distribution_upper`),
    `Image (All)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_all_mean, pred_stage_img_all_lower, pred_stage_img_all_upper),
    `2WW (All)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_all_mean, pred_stage_2ww_all_lower, pred_stage_2ww_all_upper),
    `Image (50%)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_haf_mean, pred_stage_img_haf_lower, pred_stage_img_haf_upper),
    `2WW (50%)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_haf_mean, pred_stage_2ww_haf_lower, pred_stage_2ww_haf_upper),
    `Image (25%) + 2WW (25%)`  = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_2ww_qua_mean, pred_stage_img_2ww_qua_lower, pred_stage_img_2ww_qua_upper)
  ) %>%
  print(width = Inf)

write.csv(oeso_stage, file.path(output, "MV_stagedist_oeso_stage_2026-09.csv"))


stom_stage <- stage_summary %>%
  filter(Site == "stom") %>%
  transmute(
    Site, Stage,
    `Stage distribution`      = sprintf("%.1f (%.1f-%.1f)", `Stage distribution_mean`, `Stage distribution_lower`, `Stage distribution_upper`),
    `Image (All)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_all_mean, pred_stage_img_all_lower, pred_stage_img_all_upper),
    `2WW (All)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_all_mean, pred_stage_2ww_all_lower, pred_stage_2ww_all_upper),
    `Image (50%)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_haf_mean, pred_stage_img_haf_lower, pred_stage_img_haf_upper),
    `2WW (50%)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_haf_mean, pred_stage_2ww_haf_lower, pred_stage_2ww_haf_upper),
    `Image (25%) + 2WW (25%)`  = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_2ww_qua_mean, pred_stage_img_2ww_qua_lower, pred_stage_img_2ww_qua_upper)
  ) %>%
  print(width = Inf)

write.csv(stom_stage, file.path(output, "MV_stagedist_stom_stage_2026-09.csv"))


galb_stage <- stage_summary %>%
  filter(Site == "galb") %>%
  transmute(
    Site, Stage,
    `Stage distribution`      = sprintf("%.1f (%.1f-%.1f)", `Stage distribution_mean`, `Stage distribution_lower`, `Stage distribution_upper`),
    `Image (All)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_all_mean, pred_stage_img_all_lower, pred_stage_img_all_upper),
    `2WW (All)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_all_mean, pred_stage_2ww_all_lower, pred_stage_2ww_all_upper),
    `Image (50%)`              = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_haf_mean, pred_stage_img_haf_lower, pred_stage_img_haf_upper),
    `2WW (50%)`                = sprintf("%.1f (%.1f-%.1f)", pred_stage_2ww_haf_mean, pred_stage_2ww_haf_lower, pred_stage_2ww_haf_upper),
    `Image (25%) + 2WW (25%)`  = sprintf("%.1f (%.1f-%.1f)", pred_stage_img_2ww_qua_mean, pred_stage_img_2ww_qua_lower, pred_stage_img_2ww_qua_upper)
  ) %>%
  print(width = Inf)


write.csv(galb_stage, file.path(output, "MV_stagedist_galb_stage_2026-09.csv"))

surv_cols <- c("Five-year survival", grep("^pred_surv_", colnames(psa_master), value = TRUE))

survival_summary <- psa_master %>%
  distinct(psa_run, Site, across(all_of(surv_cols))) %>%
  group_by(Site) %>%
  summarise(across(all_of(surv_cols),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        lower = ~ quantile(.x, probs = 0.025, na.rm = TRUE),
                        upper = ~ quantile(.x, probs = 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),.groups = "drop")

five_year <- survival_summary %>%
  transmute(
    Site,
    `Five-year survival`       = sprintf("%.1f (%.1f-%.1f)", `Five-year survival_mean`, `Five-year survival_lower`, `Five-year survival_upper`),
    `Image (All)`              = sprintf("%.1f (%.1f-%.1f)", pred_surv_img_all_mean, pred_surv_img_all_lower, pred_surv_img_all_upper),
    `2WW (All)`                = sprintf("%.1f (%.1f-%.1f)", pred_surv_2ww_all_mean, pred_surv_2ww_all_lower, pred_surv_2ww_all_upper),
    `Image (50%)`              = sprintf("%.1f (%.1f-%.1f)", pred_surv_img_haf_mean, pred_surv_img_haf_lower, pred_surv_img_haf_upper),
    `2WW (50%)`                = sprintf("%.1f (%.1f-%.1f)", pred_surv_2ww_haf_mean, pred_surv_2ww_haf_lower, pred_surv_2ww_haf_upper),
    `Image (25%) + 2WW (25%)`  = sprintf("%.1f (%.1f-%.1f)", pred_surv_img_2ww_qua_mean, pred_surv_img_2ww_qua_lower, pred_surv_img_2ww_qua_upper)
  ) %>%
  print(width = Inf)

write.csv(five_year, file.path(output, "MV_five_year_2026-09.csv"))