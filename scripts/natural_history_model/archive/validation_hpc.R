# Load data from HPC job script 
# Validation for parameters plus age, sex and ethnicity

rm(list = ls())
library(tidyverse)
library(scales)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))
source(file.path(scr, "natural_history_model", "fn_tp_recalibration.R"))

d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

oeso <- readRDS(file.path(output, "optim", "optim_plus_param_RS_100K_oeso.rds"))
stom <- readRDS(file.path(output, "optim", "optim_plus_param_RS_100K_stom.rds"))
panc <- readRDS(file.path(output, "optim", "optim_plus_param_RS_100K_panc.rds"))
galb <- readRDS(file.path(output, "optim", "optim_plus_param_RS_100K_galb.rds"))

# oeso <- readRDS(file.path(output, "optim_pars_RS1e+05_plus_param_oeso.rds"))
# stom <- readRDS(file.path(output, "optim_pars_RS1e+05_plus_param_stom.rds"))
# panc <- readRDS(file.path(output, "optim_pars_RS1e+05_plus_param_panc.rds"))
# galb <- readRDS(file.path(output, "optim_pars_RS1e+05_plus_param_galb.rds"))

# Load optimised parameters
# params_file0 <- read_csv(file.path(output, "optim_pars_RS100000_age_sex_ethn_HR.csv"))

params_file <- rbind(galb$final_params, oeso$final_params, panc$final_params,
                      stom$final_params) 
colnames(params_file) <- c("tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4", 
                           "beta_age70plus", "beta_female", "beta_nonwhite")
params_file <- as.data.frame(params_file)
site <- c("galb", "oeso", "panc", "stom")
params_file <- cbind(site, params_file)

# calculate distribution without counting undiagnosed cases within month limit
dist <- predict_stage_distribution(d4s, params_file, n_sim = 200, months=24,
                                   plus_param = TRUE)

# plot whole group
plot_comparison(dist)


# plot subgroups
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

ggsave(file.path(output, "subgroup_valiation3.pdf"), fig, width = 40, height = 60, units = "cm", dpi = 300)


# Present parameters
all_params <- params_file
all_params[, grepl("tp|beta", colnames(all_params))] <- round(all_params[, grepl("tp|beta", colnames(all_params))], 3)




