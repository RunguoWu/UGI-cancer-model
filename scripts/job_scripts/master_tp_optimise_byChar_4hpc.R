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

d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) %>% 
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

###### input parameters
# locate to the project root directory
args <- commandArgs(TRUE)
i <- as.numeric(args[1])
output_folder <- as.character(args[2])
char_name <- as.character(args[3])
optim_method <- as.character(args[4])
n_random_samples <- as.numeric(args[5])

print(args)
n_cores <- as.numeric(Sys.getenv('NSLOTS'))
######
combinations <- d4s %>%
  distinct(site, .data[[char_name]]) %>%
  arrange(.data[[char_name]], site)

site_val <- combinations$site[i]
char_val <- combinations[[char_name]][i]
dw_list <- dw_ugi[[site_val]]

print(paste("Processing:", "site =", site_val, ";", char_name, "=", char_val))
print(paste("Method:", optim_method, "; Sampling no.:", n_random_samples))

# Filter data for this combination
subset_data <- d4s %>%
  filter(site == site_val, .data[[char_name]] == char_val)

result <- estimate_parameters(
  patient_data = subset_data,
  n_sim = 10000, # repetitions for predicting initial distribution
  n_core = n_cores,
  seed_n = 1234,
  calib_use_month = TRUE, # likelihood consider stage X at month Y, or stage X only
  dw_list,
  plus_param = FALSE, # TRUE = plus age sex and ethn
  optim_method = optim_method,
  n_param_try = n_random_samples
)

saveRDS(result, file.path(output_folder, paste0("optim_", optim_method,
                                                "_", n_random_samples/1000, "K",  
                                                "_", char_name, char_val,  
                                                "_", site_val, ".rds")))

