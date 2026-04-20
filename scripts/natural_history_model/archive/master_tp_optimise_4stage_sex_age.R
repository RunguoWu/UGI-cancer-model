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


# Optimisation on cancer patients by site ---------------------------------
# Get unique combinations
combinations <- d4s %>%
  distinct(site) %>%
  arrange(site)

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
  beta_age70plus = numeric(),
  beta_female = numeric(),
  beta_nonwhite = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)
# saveRDS(results2, file.path(wd, "optim_pars_RS_interim_age_sex_ethn.rds"))
# Process each combination
for (i in 1:nrow(combinations)) {

  results2 <- readRDS(file.path(wd, "optim_pars_RS_interim_age_sex_ethn.rds"))
  
  site_val <- combinations$site[i]
  
  dw_list <- dw_ugi[[site_val]]
  
  print(paste("Processing:", site_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val)

  # Run optimization
  result1 <- estimate_parameters(
    patient_data = subset_data,
    n_sim = 5000, # repetitions for predicting initial distribution
    n_core = 1,
    seed_n = 1234,
    calib_use_month = TRUE, # likelihood consider stage X at month Y, or stage X only
    dw_list,
    plus_param = TRUE, # TRUE = plus age sex and ethn
    # Only work when plus_param is true
    # 80% sample from Gaussian distribution
    suggest_logit_mean = c(-0.16, -0.19, -0.09), 
    # above are means for Gaussian sampling for age70plus, female, and nonwhite,
    # converted from coefficients of AFT model on diagnostic interval
    suggest_logit_sd = c(0.1, 0.1, 0.1),
    # 20% sample from uniform distribution 
    # Boundries for uniform distribution
    lower = -1.5, 
    upper = 0.2, 
    sigma = 5, # penalise 281 for three extreme values,i.e. -1
    optim_method = "LBB",
    n_param_try = 10
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
    beta_age70plus = result2$final_params[8],
    beta_female = result2$final_params[9],
    beta_nonwhite = result2$final_params[10],
    likelihood = result2$final_likelihood,
    convergence = result2$convergence_code,
    stringsAsFactors = FALSE
  ))
  
  saveRDS(results2, file.path(wd, "optim_pars_RS_interim_age_sex_ethn.rds"))
}
# Save results
write.csv(results2, file.path(output, "optim_pars_RS100000_age_sex_ethn_HR.csv"), row.names = FALSE)


# Validation --------------------------------------------------------------

# Load optimised parameters
params_file <- read_csv(file.path(output, "optim_pars_RS100000_age_sex_ethn_HR.csv"))

# calculate distribution without counting undiagnosed cases within month limit
dist <- predict_stage_distribution(d4s, params_file, n_sim = 200, months=24,
                                   plus_param = TRUE)

plot_comparison(dist)

compare_distributions(xxx)
plot_comparison(xxx)

plot_prep_sub

plot_prep_sub(dist, sub_name)

library(gridExtra)  # or library(patchwork) if you prefer
library(grid)

# Generate the 6 plots with titles but without legends
plot_o7 <- plot_comparison(plot_prep_sub(dist, "o7")) + ggtitle("Over 70") + theme(legend.position = "none")
plot_u7 <- plot_comparison(plot_prep_sub(dist, "u7")) + ggtitle("Under 70") + theme(legend.position = "none")
plot_fe <- plot_comparison(plot_prep_sub(dist, "fe")) + ggtitle("Female") + theme(legend.position = "none")
plot_ma <- plot_comparison(plot_prep_sub(dist, "ma")) + ggtitle("Male") + theme(legend.position = "none")
plot_nw <- plot_comparison(plot_prep_sub(dist, "nw")) + ggtitle("Non-White") + theme(legend.position = "none")
plot_wh <- plot_comparison(plot_prep_sub(dist, "wh")) + ggtitle("White") + theme(legend.position = "none")

# Create one plot with legend to extract it
plot_with_legend <- plot_comparison(plot_prep_sub(dist, "o7")) + 
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Extract the legend
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(plot_with_legend)

# Combine plots and legend
fig <- grid.arrange(plot_o7, plot_u7, 
             plot_fe, plot_ma, 
             plot_nw, plot_wh, 
             legend,
             ncol = 2, nrow = 4,
             heights = c(1, 1, 1, 0.2))

ggsave(file.path(output, "subgroup_valiation.pdf"), fig, width = 30, height = 60, units = "cm", dpi = 300)

