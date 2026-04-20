# Master file for optimising parameters of the short-term natural history model

rm(list = ls())
library(tidyverse)
library(scales)
library(foreach)
library(doSNOW)
library(parallel)
library(optimParallel)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))
source(file.path(scr, "natural_history_model", "fn_tp_recalibration.R"))

dt <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds"))
d4s <- dt %>% mutate(
  age_group3 = case_when(age_index < 60 ~ "<60",
                         age_index >= 60 & age_index < 80 ~ "60-79",
                         age_index >=80 ~ ">=80"
  ),
  age_group3 = factor(age_group3, levels = c("<60", "60-79", ">=80")),
  month = ceiling(time2diag/30.5),
  diagnosed_stage = as.integer(stage_imp)
) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age_group3)


# Optimisation on cancer patients by site ---------------------------------
# Get unique combinations
combinations <- d4s %>%
  distinct(site) %>%
  arrange(site)

############################
## L-BFGS-B method ----
############################
# Initialize results dataframe
results1 <- data.frame(
  site = character(),
  n_patients = integer(),
  tp12 = numeric(),
  tp23 = numeric(),
  tp34 = numeric(),
  tp1 = numeric(),
  tp2 = numeric(),
  tp3 = numeric(),
  tp4 = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)

saveRDS(results1, file.path(wd, "optim_pars_BFGS_interim.rds"))

# Process each combination
for (i in 1:nrow(combinations)) {
  
  results1 <- readRDS(file.path(wd, "optim_pars_BFGS_interim.rds"))
  
  site_val <- combinations$site[i]
  
  print(paste("Processing:", site_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val)

  # Run optimization
  result1 <- estimate_parameters_BFGS(
    patient_data = subset_data,
    lower_params = bounds$lower,
    upper_params = bounds$upper,
    # save_interval = NA,
    filename_base = paste0("opt_", site_val),
    save_final = FALSE,
    n_initial_params = 50,
    n_core = 8
  )

  # Store results
  results1 <- rbind(results1, data.frame(
    site = site_val,
    n_patients = nrow(subset_data),
    tp12 = result1$final_params[1],
    tp23 = result1$final_params[2],
    tp34 = result1$final_params[3],
    tp1 = result1$final_params[4],
    tp2 = result1$final_params[5],
    tp3 = result1$final_params[6],
    tp4 = result1$final_params[7],
    likelihood = result1$final_likelihood,
    convergence = result1$convergence_code,
    stringsAsFactors = FALSE
  ))
  
  saveRDS(results1, file.path(wd, "optim_pars_BFGS_interim.rds"))
  
}
# Save results
write.csv(results1, file.path(output, "optim_pars_BFGS.csv"), row.names = FALSE)

############################
## Random Search method ----
############################
results2 <- data.frame(
  site = character(),
  n_patients = integer(),
  tp12 = numeric(),
  tp23 = numeric(),
  tp34 = numeric(),
  tp1 = numeric(),
  tp2 = numeric(),
  tp3 = numeric(),
  tp4 = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)
# saveRDS(results2, file.path(wd, "optim_pars_RS_interim.rds"))
# Process each combination
for (i in 1:nrow(combinations)) {
  
  results2 <- readRDS(file.path(wd, "optim_pars_RS_interim.rds"))
  
  site_val <- combinations$site[i]
  
  print(paste("Processing:", site_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val)
  
  # Run optimization
  result2 <- estimate_parameters_rs(
    patient_data = subset_data,
    lower_params = bounds$lower,
    upper_params = bounds$upper,
    save_interval = NA,
    filename_base = paste0("opt_", site_val),
    save_final = FALSE, 
    n_random_samples = 100000,
    n_core = 8
  )
  
  # Store results
  results2 <- rbind(results2, data.frame(
    site = site_val,
    n_patients = nrow(subset_data),
    tp12 = result2$final_params[1],
    tp23 = result2$final_params[2],
    tp34 = result2$final_params[3],
    tp1 = result2$final_params[4],
    tp2 = result2$final_params[5],
    tp3 = result2$final_params[6],
    tp4 = result2$final_params[7],
    likelihood = result2$final_likelihood,
    convergence = result2$convergence_code,
    stringsAsFactors = FALSE
  ))
  
  saveRDS(results2, file.path(wd, "optim_pars_RS_interim.rds"))
}
# Save results
write.csv(results2, file.path(output, "optim_pars_RS.csv"), row.names = FALSE)

############################
## Nelder-Mead method----
############################
results3 <- data.frame(
  site = character(),
  n_patients = integer(),
  tp12 = numeric(),
  tp23 = numeric(),
  tp34 = numeric(),
  tp1 = numeric(),
  tp2 = numeric(),
  tp3 = numeric(),
  tp4 = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)
# saveRDS(results3, file.path(wd, "optim_pars_NM_interim.rds"))
# Process each combination
for (i in 1:nrow(combinations)) {
  
  results3 <- readRDS(file.path(wd, "optim_pars_NM_interim.rds"))
  
  site_val <- combinations$site[i]
  
  print(paste("Processing:", site_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val)
  
  result3 <- estimate_parameters_NM(
    patient_data = subset_data,
    filename_base = paste0("opt_", site_val),
    save_final = FALSE, 
    n_initial_params = 50,
    n_core = 8
  )
  
  results3 <- rbind(results3, data.frame(
    site = site_val,
    n_patients = nrow(subset_data),
    tp12 = result3$final_params[1],
    tp23 = result3$final_params[2],
    tp34 = result3$final_params[3],
    tp1 = result3$final_params[4],
    tp2 = result3$final_params[5],
    tp3 = result3$final_params[6],
    tp4 = result3$final_params[7],
    likelihood = result3$final_likelihood,
    convergence = result3$convergence_code,
    stringsAsFactors = FALSE
  ))
  
  saveRDS(results3, file.path(wd, "optim_pars_NM_interim.rds"))
}
# Save results
write.csv(results3, file.path(output, "optim_pars_NM.csv"), row.names = FALSE)



# Optimisation on cancer patients by site and age -------------------------

# Set optimization bounds
bounds <- get_parameter_bounds()

# Initialize results dataframe
results <- data.frame(
  site = character(),
  age_group3 = character(),
  n_patients = integer(),
  tp12 = numeric(),
  tp23 = numeric(),
  tp34 = numeric(),
  tp1 = numeric(),
  tp2 = numeric(),
  tp3 = numeric(),
  tp4 = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)

# Get unique combinations
combinations <- d4s %>%
  distinct(site, age_group3) %>%
  arrange(site, age_group3)

print(paste("Processing", nrow(combinations), "combinations"))

# Process each combination
for (i in 1:nrow(combinations)) {
  site_val <- combinations$site[i]
  age_val <- combinations$age_group3[i]
  
  print(paste("Processing:", site_val, "-", age_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val, age_group3 == age_val)
  
  # Skip if too few patients
  if (nrow(subset_data) < 10) {
    print(paste("Skipping - only", nrow(subset_data), "patients"))
    next
  }
  
  # Get suggested starting parameters
  initial_params <- suggest_starting_parameters(subset_data)
  
  # Run optimization
  result <- estimate_parameters(
    patient_data = subset_data,
    initial_params = initial_params,
    lower_params = bounds$lower,
    upper_params = bounds$upper,
    save_interval = NA,
    filename_base = paste0("opt_", site_val, "_", age_val),
    save_final = FALSE
  )
  
  # Store results
  results <- rbind(results, data.frame(
    site = site_val,
    age_group3 = age_val,
    n_patients = nrow(subset_data),
    tp12 = result$final_params[1],
    tp23 = result$final_params[2],
    tp34 = result$final_params[3],
    tp1 = result$final_params[4],
    tp2 = result$final_params[5],
    tp3 = result$final_params[6],
    tp4 = result$final_params[7],
    likelihood = result$final_likelihood,
    convergence = result$convergence_code,
    stringsAsFactors = FALSE
  ))
}

# Display results
print("Optimization Results:")
print(results)

# Save results
write.csv(results, file.path(output, "optimization_results_4stage_bySiteAge.csv"), row.names = FALSE)


# Validation --------------------------------------------------------------

# Load optimised parameters
# params <- read_csv(file.path(output, "optimization_results_4stage_bySite.csv"))
# params <- read_csv(file.path(output, "optimization_results_4stage_bySite_noRedFlag.csv"))

# params_file <- read_csv(file.path(output, "optim_pars_NM.csv"))
params_file <- read_csv(file.path(output, "optim_pars_RS.csv"))

# calculate distribution without counting undiagnosed cases within month limit
dist <- predict_stage_distribution(d4s, params_file, n_sim = 100, months=24)

plot_comparison(dist)


# Recalibration using total likelihood ------------------------------------
# The result was not good
# Stop using it
# 
# results_recali <- data.frame(
#   site = character(),
#   n_patients = integer(),
#   tp12 = numeric(),
#   tp23 = numeric(),
#   tp34 = numeric(),
#   tp1 = numeric(),
#   tp2 = numeric(),
#   tp3 = numeric(),
#   tp4 = numeric(),
#   likelihood = numeric(),
#   convergence = integer(),
#   stringsAsFactors = FALSE
# )
# 
# params_file <- read_csv(file.path(output, "optim_pars_NM.csv"))
# 
# for (st in names(dist)) {
#   
#   patient_data <- d4s[d4s$site == st, ]
#   
#   optimized_params <- as.numeric(params_file[params_file$site==st, grepl("tp", colnames(params_file))])
#   
#   rt <- estimate_parameters_recali(patient_data,
#                                    optimized_params,
#                                    filename_base = "mle_results_recali",
#                                    save_final = FALSE, 
#                                    seed_n =1234
#   ) 
#   
#   results_recali <- rbind(results_recali, data.frame(
#     site = st,
#     n_patients = nrow(patient_data),
#     tp12 = rt$final_params[1],
#     tp23 = rt$final_params[2],
#     tp34 = rt$final_params[3],
#     tp1 = rt$final_params[4],
#     tp2 = rt$final_params[5],
#     tp3 = rt$final_params[6],
#     tp4 = rt$final_params[7],
#     likelihood = rt$final_likelihood,
#     convergence = rt$convergence_code,
#     stringsAsFactors = FALSE
#   ))
# 
# }
# # Save results
# write.csv(results_recali, file.path(output, "optim_pars_recali.csv"), row.names = FALSE)
# 
# dist2 <- predict_stage_distribution(d4s, results_recali, n_sim = 500)
# 
# plot_comparison(dist2)


# Direct Re-Calibration ---------------------------------------------------
# Load originally optimised parameters
# params <- read_csv(file.path(output, "optimization_results_4stage_bySite.csv"))

# params_file <- read_csv(file.path(output, "optim_pars_NM.csv"))
params_file <- read_csv(file.path(output, "optim_pars_RS.csv"))
dist <- predict_stage_distribution2(d4s, params_file, n_sim = 100, months = 24)
compare_distributions(dist)
plot_comparison(dist)
# Initial prediction of distribution on initial param
# checked n_sim = 100 is enough, no need for 500

calibration_results <- iterative_calibration(
  patient_data = d4s,
  initial_params = params_file,
  max_iterations = 10,
  convergence_threshold = 0.025,
  adjustment_weight = 0.3,
  n_sim = 100,
  adaptive_weight = TRUE
)

# saveRDS(calibration_results, file.path(wd, "calibration_results_NM.rds"))
saveRDS(calibration_results, file.path(wd, "calibration_results_RS.rds"))

calibration_results <- readRDS(file.path(wd, "calibration_results_RS.rds"))

# Plot
plot_comparison(calibration_results$final_comparison, re_cali = T)


# finally check likelihood
st <- "panc"
dtt <- d4s[d4s$site==st, ]

optimized_params <- as.numeric(params_file[params_file$site==st, grepl("tp", colnames(params_file))])

recali_params <- calibration_results$calibrated_params

final_params <- as.numeric(recali_params[recali_params$site==st, grepl("tp", colnames(recali_params))])

calculate_cohort_likelihood(optimized_params, dtt, use_month = T) 
calculate_cohort_likelihood(final_params, dtt, use_month = T) 


# Intervention effect -----------------------------------------------------

recali_file <- readRDS(file.path(wd, "calibration_results_RS.rds"))

params <- recali_file$calibrated_params

pred_dist <- predict_stage_distribution2(d4s, params, 
                                         n_sim = 100, months=24, 
                                         tx=1, tx_prob=1
)

pred_dist_interv <- predict_stage_distribution2(d4s, params, 
                                                 n_sim = 100, months=24, 
                                                 tx=1.8, tx_prob=0.5
)

comparison_df <- compare_distributions_interv(pred_dist, pred_dist_interv)

plot_comparison_interv(comparison_df)



