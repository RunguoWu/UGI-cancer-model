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

# d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) %>% 
d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_upd2026.rds")) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
# opt_tag <- "opt20251114_RS_500K"
# opt_tag <- "opt20251202_RS_500K"
opt_tag <- "opt_upd2026_RS_500K"

# Fig1: Internal validation -----------------------------------------------
dist_list <- readRDS(file.path(wd, "predicted_distribution_by2Char_upd2026.rds"))

plot_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    sub_dist <- dist_list[[sub_name]] 
    
    title1 <- if (char_value1 == 0) "Under 70" else "Over 70"
    title2 <- if (char_value2 == 0) "Men" else "Women"
    title12 <- paste(title1,";", title2)
    
    sub_plot <- plot_comparison(sub_dist) + ggtitle(title12) + 
      theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
    
    plot_list[[sub_name]] <- sub_plot
  }
}

# Create one plot with legend to extract it
plot_with_legend <- plot_comparison(dist_list[[1]]) + 
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

# ggsave(file.path(output, paste0(opt_tag, "_valiation_by2Char_20251202_alls4.png")), 
#        fig, width = 30, height = 30, units = "cm", dpi = 300)
ggsave(file.path(output, paste0(opt_tag, "_valiation_by2Char_upd2026_alls4.png")), 
       fig, width = 30, height = 30, units = "cm", dpi = 300)


# Fig 2: Imaging scenario -------------------------------------------------
tx_name <- "imaging"
# tx_name <- "2ww"
# ind_output_list <- readRDS(file.path(output, paste0("ind_interv_", tx_name, "_by2Char_upd2026.rds")))
ind_output_list <- readRDS(file.path(output, paste0("ind_interv_", tx_name, "_by2Char_upd2026_bySiteOnly.rds")))

# Make plot
plot_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    comparison_df <- compare_distributions_interv(ind_output_list[[sub_name]])
    
    title1 <- if (char_value1 == 0) "Under 70" else "Over 70"
    title2 <- if (char_value2 == 0) "Men" else "Women"
    title12 <- paste(title1,";", title2)
    
    sub_plot <- plot_comparison_interv(comparison_df) + ggtitle(title12) + 
      theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
    
    plot_list[[sub_name]] <- sub_plot
  }
}


# Create one plot with legend to extract it
comparison_df <- compare_distributions_interv(ind_output_list[[1]])

plot_with_legend <- plot_comparison_interv(comparison_df) + 
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

# ggsave(file.path(output, paste0("interv_", tx_name, "_by2Char_upd2026.png")), fig, width = 28, height = 35, units = "cm", dpi = 300)
ggsave(file.path(output, paste0("interv_", tx_name, "_by2Char_upd2026_bySiteOnly.png")), fig, width = 28, height = 35, units = "cm", dpi = 300)

