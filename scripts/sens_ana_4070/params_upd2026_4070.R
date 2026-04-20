# Load data from HPC
# Validation for by-age and by-sex parameters

rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))
source(file.path(scr, "natural_history_model", "fn_tp_recalibration.R"))

d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_upd2026.rds")) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  filter(age_index >=40) %>% 
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
opt_tag <- "opt_upd2026_sa4070_RS_500K"

dt_list <- list()
for (site in site_name) {
  
  dt_list[[paste0(site, "_", char_name[1], 0, "_", char_name[2], 0)]] <- 
    readRDS(file.path(output, "optim", 
                      paste0(opt_tag, "_", 
                             char_name[1], 0, "_", 
                             char_name[2], 0, "_", site, ".rds")))
  
  dt_list[[paste0(site, "_", char_name[1], 0, "_", char_name[2], 1)]] <- 
    readRDS(file.path(output, "optim", 
                      paste0(opt_tag, "_",  
                             char_name[1], 0, "_", 
                             char_name[2], 1, "_", site, ".rds")))
  
  dt_list[[paste0(site, "_", char_name[1], 1, "_", char_name[2], 0)]] <- 
    readRDS(file.path(output, "optim", 
                      paste0(opt_tag, "_", 
                             char_name[1], 1, "_", 
                             char_name[2], 0, "_", site, ".rds")))
  
  dt_list[[paste0(site, "_", char_name[1], 1, "_", char_name[2], 1)]] <- 
    readRDS(file.path(output, "optim", 
                      paste0(opt_tag, "_", 
                             char_name[1], 1, "_", 
                             char_name[2], 1, "_", site, ".rds")))
}


## New method ----
# Use the median of the distribution of the parameters with likelihood above the threshold
params_mat <- matrix(NA, nrow = 16, ncol=7)
df_list <- list()
rownames(params_mat) <- names(dt_list)
for (sub_name in names(dt_list)) {
  records <- dt_list[[sub_name]]$params_record
  llh_records <- sapply(records, "[[", "value_cur")
  par_records <- sapply(records, "[[", "par_cur")
  
  params <- dt_list[[sub_name]]$final_params
  
  threshold <- quantile(llh_records, probs = 0.01) # .01 corresponding to 5000, with 500K iterations
  keep_indices <- which(llh_records <= threshold)
  par_records_keep <- par_records[, keep_indices]
  
  params_med <- as.numeric(apply(par_records_keep, 1, quantile, probs = 0.5))
  params_95ci <- t(apply(par_records_keep, 1, quantile, probs = c(0.025, 0.975)))
  # params_med <- as.numeric(apply(par_records_keep, 1, mean))
  
  params_mat[rownames(params_mat) == sub_name, ] <- params_med
  
  df <- data.frame(
    index = c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4"),
    estimate = params_med, # params,
    lower = params_95ci[, "2.5%"],
    upper = params_95ci[, "97.5%"]
  )
  df_list[[sub_name]] <- df
}

params_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    params <- params_mat[grepl(sub_name, rownames(params_mat)),  ]
    
    site <- sub("_.*", "", rownames(params))
    params <- as.data.frame(params)
    params <- cbind(site, params)
    colnames(params) <- c("site", "tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
    rownames(params) <- NULL
    
    params_list[[sub_name]] <- params
  }
}

# Present parameters
all_params <- do.call(rbind, params_list)
all_params[, grepl("tp", colnames(all_params))] <- round(all_params[, grepl("tp", colnames(all_params))], 3)
subgroups <- c(rep("<70; Men", 4), rep("<70; Women", 4), rep(">=70; Men", 4), rep(">=70; Women", 4))
all_params <- cbind(subgroups, all_params)
rownames(all_params) <- NULL
x <- all_params %>% arrange(site, subgroups) %>% relocate(site)
write.csv(x, file.path(output, "ugi-cpdi model parameters_upd2026_sa4070.csv"))