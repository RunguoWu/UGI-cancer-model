rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))

d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
opt_tag <- "opt20251202_RS_500K"


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

plot_list <- list()
df_list <- list()
for (subgroup in names(dt_list)) {
  records <- dt_list[[subgroup]]$params_record
  llh_records <- sapply(records, "[[", "value_cur")
  par_records <- sapply(records, "[[", "par_cur")
  
  params <- dt_list[[subgroup]]$final_params
  
  threshold <- quantile(llh_records, probs = 0.01)
  keep_indices <- which(llh_records <= threshold)
  par_records_keep <- par_records[, keep_indices]
  params_95ci <- t(apply(par_records_keep, 1, quantile, probs = c(0.025, 0.975)))
  params_med <- as.numeric(apply(par_records_keep, 1, quantile, probs = 0.5))
  
  df <- data.frame(
    index = 1:7,
    estimate = params_med, # params,
    lower = params_95ci[, "2.5%"],
    upper = params_95ci[, "97.5%"]
  )
  
  fig <- ggplot(df, aes(x = index, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "Parameter Index", 
         y = "Value",
         title = "Parameters with 95% Confidence Intervals") +
    theme_minimal()
  
  df_list[[subgroup]] <- df
  plot_list[[subgroup]] <- fig
}

saveRDS(df_list, file.path(wd, "params_list_20251202_withCI.rds"))














