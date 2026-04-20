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
d4s <- dt %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

# Dwelling years
# Cite https://doi.org/10.1371/journal.pone.0279227
# all stage 4 use 0.25-1
oeso_dw <- c(stage1_lower = 2, stage1_upper = 5,
             stage2_lower = 0.5, stage2_upper = 2,
             stage3_lower = 0.5, stage3_upper = 2,
             stage4_lower = 0.25, stage4_upper =1)

stom_dw <- c(stage1_lower = 2, stage1_upper = 5,
             stage2_lower = 0.5, stage2_upper = 2,
             stage3_lower = 0.5, stage3_upper = 2,
             stage4_lower = 0.25, stage4_upper =1)

# in literature, some ranges = <1-<1, assume 0.25-1
galb_dw <- c(stage1_lower = 0.5, stage1_upper = 3,
             stage2_lower = 0.25, stage2_upper = 1,
             stage3_lower = 0.25, stage3_upper = 1,
             stage4_lower = 0.25, stage4_upper =1)

panc_dw <- c(stage1_lower = 0.5, stage1_upper = 2,
             stage2_lower = 0.5, stage2_upper = 2,
             stage3_lower = 0.25, stage3_upper = 1,
             stage4_lower = 0.25, stage4_upper =1)

dw_ugi <- list("oeso"=oeso_dw, "stom"=stom_dw, "galb"=galb_dw, "panc"=panc_dw)

char_name <- "female"

# Optimisation on cancer patients by age70plus and site -------------------
# Get unique combinations
combinations <- d4s %>%
  distinct(site, .data[[char_name]]) %>%
  arrange(.data[[char_name]], site)

results <- data.frame(
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

for (i in 1:nrow(combinations)) {
  
  # results <- readRDS(file.path(wd, "optim_pars_RS_interim_age_sex_ethn.rds"))
  
  site_val <- combinations$site[i]
  char_val <- combinations[[char_name]][i]
  
  dw_list <- dw_ugi[[site_val]]
  
  print(paste("Processing:", "site =", site_val, ";", char_name, "=", char_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val, .data[[char_name]] == char_val)
  
  # Run optimization
  
  # patient_data = subset_data
  # n_sim = 5000
  # n_core = 2
  # seed_n = 1234
  # calib_use_month = TRUE
  # plus_param = TRUE
  # suggest_logit_mean = c(-0.16, -0.19, -0.09) 
  # suggest_logit_sd = c(0.1, 0.1, 0.1)
  # lower = -1.5 
  # upper = 0.2
  # sigma = 5
  # optim_method = "NM"
  # n_param_try = 10
  
  result <- estimate_parameters(
    patient_data = subset_data,
    n_sim = 10000, # repetitions for predicting initial distribution
    n_core = 2,
    seed_n = 1234,
    calib_use_month = TRUE, # likelihood consider stage X at month Y, or stage X only
    dw_list,
    plus_param = FALSE, # TRUE = plus age sex and ethn
    optim_method = "RS",
    n_param_try = 1000
  )
  
  # Store results
  results <- rbind(results, data.frame(
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
  
  saveRDS(results, file.path(wd, "optim_pars_RS_interim_age_sex_ethn.rds"))
}
# Save results
write.csv(results, file.path(output, "optim_pars_RS100000_age_sex_ethn_HR.csv"), row.names = FALSE)
