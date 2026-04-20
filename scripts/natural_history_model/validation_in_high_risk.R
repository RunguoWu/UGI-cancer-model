# Validate the treatment effect of imaging and urgent referral recommendation

rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)  # or library(patchwork) if you prefer
library(grid)
library(survival)
library(rstpm2)
library(zoo)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))
source(file.path(scr, "natural_history_model", "fn_tp_recalibration.R"))

# Patient data
# dt <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds"))
# params_list <- readRDS(file.path(wd, "params_list_20251202.rds"))
# tx_list <- readRDS(file.path(wd, "tx_list.rds"))

dt <- readRDS(file.path(wd, "study_pop_symp_stageImputed_upd2026.rds")) 
dt_na <- dt %>% 
  filter(ng12_red_flag == "No red flag") 

params_list <- readRDS(file.path(wd, "params_list_upd2026.rds"))

# Load tx
# tx_list <- readRDS(file.path(wd, "tx_list_upd2026.rds"))
tx_list <- readRDS(file.path(wd, "tx_list_upd2026_bySiteOnly.rds"))

# Subgroups
site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
tx_name <- "2ww" # "imaging" #     

if (tx_name=="2ww") d4s <- dt %>% filter(ng12_red_flag2 == "2 Week Wait")
if (tx_name=="imaging") d4s <- dt %>% filter(ng12_red_flag2 == "Imaging")


ind_output_list <- list()

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    params <- params_list[[sub_name]]
    
    tx <- sapply(tx_list[[sub_name]], "[", tx_name)
    names(tx) <- gsub(paste0("\\.", tx_name, "$"), "", names(tx))
    
    # Load average index stage distribution by site and subgroups
    # index_dist <- read.csv(file.path(output, paste0("initial_dist_", sub_name, ".csv")))
    # this is too broad, use by month distribution
    
    sub_data_na <- dt_na %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    initial_stage <- list()
    dist_list <- estimate_starting_distribution(sub_data_na, 
                                                params, 
                                                plus_param=FALSE, 
                                                average_dist = FALSE)
    for (st in site_name) {
      initial_stage[[st]] <- dist_list[[st]] %>% 
        left_join(sub_data_na[,c("e_patid", "month", "diagnosed_stage")]) %>%
        group_by(diagnosed_stage, month) %>% 
        summarise(X1 = mean(X1), X2 = mean(X2), X3 = mean(X3), X4 = mean(X4)) %>% 
        mutate(site = st) %>% 
        relocate(site) %>% 
        ungroup()
    }
    initial_stage <- as.data.frame(do.call(rbind, initial_stage))
    
    # some month values are missing in initial_stage
    # impute them using neighbouring values
    complete_grid <- expand.grid(
      site = c("galb", "stom", "panc", "oeso"),
      diagnosed_stage = 1:4,
      month = 1:24
    )
    initial_stage_complete <- complete_grid %>%
      left_join(initial_stage, by = c("site", "diagnosed_stage", "month"))
    
    initial_stage_imputed <- initial_stage_complete %>%
      arrange(site, diagnosed_stage, month) %>%
      group_by(site, diagnosed_stage) %>%
      mutate(
        X1 = zoo::na.approx(X1, na.rm = FALSE),
        X2 = zoo::na.approx(X2, na.rm = FALSE),
        X3 = zoo::na.approx(X3, na.rm = FALSE),
        X4 = zoo::na.approx(X4, na.rm = FALSE)
      ) %>%
      fill(X1, X2, X3, X4, .direction = "up") %>% # fills NAs with the next non-NA value
      fill(X1, X2, X3, X4, .direction = "down") %>% 
      ungroup()
      
    pred_dist_interv <- predict_stage_distribution2(sub_data, params, 
                                                    n_sim = 1000, months=24, 
                                                    tx=tx, tx_prob=1,
                                                    use_avg_start_dist = TRUE,
                                                    avg_start_dist = initial_stage_imputed
                                                    )
    
    ind_output_list[[sub_name]] <- pred_dist_interv
  }
}

# Make plot
plot_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    comparison_df <- compare_distributions_interv2(ind_output_list[[sub_name]])
    
    title1 <- if (char_value1 == 0) "Under 70" else "Over 70"
    title2 <- if (char_value2 == 0) "Men" else "Women"
    title12 <- paste(title1,";", title2)
    
    sub_plot <- plot_comparison_interv2(comparison_df) + ggtitle(title12) + 
      theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
    
    plot_list[[sub_name]] <- sub_plot
  }
}

# Create one plot with legend to extract it
comparison_df <- compare_distributions_interv2(ind_output_list[[1]])

plot_with_legend <- plot_comparison_interv2(comparison_df) + 
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
fig <- grid.arrange(plot_list[[1]], plot_list[[2]], 
                    plot_list[[3]], plot_list[[4]], 
                    legend,
                    ncol = 2, nrow = 3,
                    heights = c(1, 1, 0.2))

# ggsave(file.path(output, paste0("valid_high_risk_", tx_name, ".png")), fig, width = 28, height = 35, units = "cm", dpi = 300)
# ggsave(file.path(output, paste0("valid_high_risk_", tx_name, "_upd2026.png")), fig, width = 28, height = 35, units = "cm", dpi = 300)
ggsave(file.path(output, paste0("valid_high_risk_", tx_name, "_upd2026_bySiteOnly.png")), fig, width = 28, height = 35, units = "cm", dpi = 300)



