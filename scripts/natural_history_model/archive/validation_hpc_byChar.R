# Load data from HPC
# Validation for by-age and by-sex parameters

rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)  # or library(patchwork) if you prefer
library(grid)

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

site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")

dt_list <- list()

for (site in site_name) {
  for (char in char_name) {
    
    dt_list[[paste0(site, "_", char, "0")]] <- 
      readRDS(file.path(output, "optim", paste0("optim_RS_100K_", char, 0, "_", site, ".rds")))
    dt_list[[paste0(site, "_", char, "1")]] <- 
      readRDS(file.path(output, "optim", paste0("optim_RS_100K_", char, 1, "_", site, ".rds")))
  }
}

# check
library(patchwork)

# Create a list of ggplot objects
plot_list <- lapply(seq_along(dt_list), function(i) {
  df <- data.frame(subgroup = seq_along(dt_list[[i]]$value_convergence_record),
                   value = dt_list[[i]]$value_convergence_record)
  ggplot(df, aes(x = subgroup, y = value)) +
    geom_line() +
    ggtitle(paste("subgroup", i)) +
    theme_minimal(base_size = 10)
})

# Combine into a 4x4 grid
final_fig <- wrap_plots(plot_list, ncol = 4)
final_fig

ggsave(file.path(output,  "convergence_4x4.png"), final_fig, width = 12, height = 10)

sub_list <- list()
for (char in char_name) {
  
  dt_list_sub0 <- dt_list[grepl(paste0(char,0), names(dt_list))]
  dt_list_sub1 <- dt_list[grepl(paste0(char,1), names(dt_list))]
  
  params0 <- do.call(rbind,lapply(dt_list_sub0, function(x) x$final_params))
  params1 <- do.call(rbind,lapply(dt_list_sub1, function(x) x$final_params))

  site0 <- sub("_.*", "", rownames(params0))
  site1 <- sub("_.*", "", rownames(params1))
  
  params0 <- as.data.frame(params0)
  params0 <- cbind(site0, params0)
  
  params1 <- as.data.frame(params1)
  params1 <- cbind(site, params1)
  
  colnames(params0) <- c("site", "tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
  colnames(params1) <- c("site", "tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
  
  rownames(params0) <- NULL
  rownames(params1) <- NULL
  
  sub_list[[paste0(char,0)]] <- params0
  sub_list[[paste0(char,1)]] <- params1
}

dist_a0 <- predict_stage_distribution(subset(d4s, age70plus==0), 
                                      sub_list[["age70plus0"]], n_sim = 200, months=24,
                                      plus_param = FALSE)


dist_a1 <- predict_stage_distribution(subset(d4s, age70plus==1), 
                                      sub_list[["age70plus1"]], n_sim = 200, months=24,
                                      plus_param = FALSE)


dist_f0 <- predict_stage_distribution(subset(d4s, female==0), 
                                      sub_list[["female0"]], n_sim = 200, months=24,
                                      plus_param = FALSE)


dist_f1 <- predict_stage_distribution(subset(d4s, female==1), 
                                      sub_list[["female1"]], n_sim = 200, months=24,
                                      plus_param = FALSE)

plot_a0 <- plot_comparison(dist_a0) + ggtitle("Under 70") + 
  theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
plot_a1 <- plot_comparison(dist_a1) + ggtitle("Over 70") + 
  theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
plot_f0 <- plot_comparison(dist_f0) + ggtitle("Male") + 
  theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
plot_f1 <- plot_comparison(dist_f1) + ggtitle("Female") + 
  theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")

# Create one plot with legend to extract it
plot_with_legend <- plot_comparison(dist_f0) + 
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
fig <- grid.arrange(plot_a0, plot_a1, 
                    plot_f0, plot_f1, 
                    legend,
                    ncol = 2, nrow = 3,
                    heights = c(1, 1, 0.2))

ggsave(file.path(output, "valiation_byChar.png"), fig, width = 40, height = 50, units = "cm", dpi = 300)

# Present parameters
all_params <- do.call(rbind, sub_list)
all_params[, grepl("tp", colnames(all_params))] <- round(all_params[, grepl("tp", colnames(all_params))], 3)

subgroups <- c(rep("<70 yo", 4), rep(">=70 yo", 4), rep("Men", 4), rep("Women", 4))
all_params <- cbind(subgroups, all_params)
rownames(all_params) <- NULL
