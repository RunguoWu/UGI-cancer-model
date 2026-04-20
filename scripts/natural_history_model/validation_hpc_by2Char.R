# Load data from HPC
# Validation for by-age and by-sex parameters

rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)
library(writexl)

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

# check convergence -------------------------------------------------------
# Create a list of ggplot objects
library(patchwork)
plot_list <- lapply(seq_along(dt_list), function(i) {
  df <- data.frame(subgroup = seq_along(dt_list[[i]]$value_convergence_record),
                   value = dt_list[[i]]$value_convergence_record)
  ggplot(df, aes(x = subgroup, y = value)) +
    geom_line() +
    ggtitle(names(dt_list[i])) +
    theme_minimal(base_size = 10)
})

# Combine into a 4x4 grid
final_fig <- wrap_plots(plot_list, ncol = 4)
final_fig

ggsave(file.path(output, paste0(opt_tag, "convergence2Char_4x4.png")), final_fig, width = 12, height = 10)


# prepare for validation --------------------------------------------------

## Old method ----
# Use the parameters with the maximum likelihood
params_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    dt_list_sub <- dt_list[grepl(sub_name, names(dt_list))]
    
    params <- do.call(rbind,lapply(dt_list_sub, function(x) x$final_params))
    
    site <- sub("_.*", "", rownames(params))
    params <- as.data.frame(params)
    params <- cbind(site, params)
    colnames(params) <- c("site", "tp12", "tp23", "tp34", "tp1", "tp2", "tp3", "tp4")
    rownames(params) <- NULL
    
    params_list[[sub_name]] <- params
  }
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
# save parameter uncertainty 95% CI
# saveRDS(df_list, file.path(wd, "params_list_20251202_withCI.rds"))
saveRDS(df_list, file.path(wd, "params_list_upd2026_withCI.rds"))

df_list <- readRDS(file.path(wd, "params_list_upd2026_withCI.rds"))

df_expo <- do.call(rbind, df_list) %>% 
  mutate(est_ci = paste0(round(estimate, 3), " (", round(lower, 3), "-", 
                         round(upper, 3), ")")) %>% 
  select(-estimate, -lower, -upper) %>% 
  rownames_to_column("id") %>% 
  mutate(id = gsub("\\.\\d+$", "", id)) %>%
  pivot_wider(
    id_cols     = id,
    names_from  = index,
    values_from = est_ci   
  ) %>% 
  mutate(
    cancer_site = case_when(
      startsWith(id, "galb") ~ "Gallbladder",
      startsWith(id, "oeso") ~ "Oesophagus",
      startsWith(id, "stom") ~ "Stomach",
      startsWith(id, "panc") ~ "Pancreas"
    ),
    age_sex = case_when(
      grepl("age70plus0_female0", id) ~ "<70; Men",
      grepl("age70plus0_female1", id) ~ "<70; Women",
      grepl("age70plus1_female0", id) ~ "\u226570; Men",
      grepl("age70plus1_female1", id) ~ "\u226570; Women"
    )
  ) %>% 
  select(cancer_site, age_sex, everything(), -id) %>%
  arrange(factor(cancer_site, levels = c("Pancreas", "Oesophagus", "Stomach", "Gallbladder")))

write_xlsx(df_expo, file.path(output, "ugi_cpdi_param_CI_upd2026.xlsx"))

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

# saveRDS(params_list, file.path(wd, "params_list_20251202.rds"))
saveRDS(params_list, file.path(wd, "params_list_upd2026.rds"))

# Present parameters
all_params <- do.call(rbind, params_list)
all_params[, grepl("tp", colnames(all_params))] <- round(all_params[, grepl("tp", colnames(all_params))], 3)
subgroups <- c(rep("<70; Men", 4), rep("<70; Women", 4), rep(">=70; Men", 4), rep(">=70; Women", 4))
all_params <- cbind(subgroups, all_params)
rownames(all_params) <- NULL
x <- all_params %>% arrange(site, subgroups) %>% relocate(site)
write.csv(x, file.path(output, "ugi-cpdi model parameters_upd2026.csv"))

# Simulate for validation -------------------------------------------------
# Simulate initial stage distribution by site and subgroups
# params_list <- readRDS(file.path(wd, "params_list_20251202.rds"))
params_list <- readRDS(file.path(wd, "params_list_upd2026.rds"))

initial_stage <- list()
initial_stage2 <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    params <- params_list[[sub_name]]
    
    dist_list <- estimate_starting_distribution(sub_data, 
                                                params, 
                                                plus_param=FALSE, 
                                                average_dist = FALSE)
    for (st in site_name) {
      initial_stage[[sub_name]][[st]] <- dist_list[[st]] %>% 
        left_join(sub_data[,c("e_patid", "diagnosed_stage")]) %>%
        group_by(diagnosed_stage) %>% 
        summarise(mean(X1), mean(X2), mean(X3), mean(X4)) %>% 
        mutate(site = st) %>% 
        relocate(site)
    }
    initial_stage2[[sub_name]] <- do.call(rbind, initial_stage[[sub_name]])
  }
}

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)

    write.csv(initial_stage2[[sub_name]], 
              file.path(output, paste0("initial_dist_", sub_name, ".csv")))
  }
}

dist_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    params <- params_list[[sub_name]]
    
    # scale up tp to see results
    # this works well. 
    # simpler than recalibration method below.
    # params[, c("tp12", "tp23", "tp34")] <- params[, c("tp12", "tp23", "tp34")]*1.2

    
    sub_dist <- predict_stage_distribution(sub_data,
                                           params, n_sim = 200, months=24,
                                           plus_param = FALSE)
    
    dist_list[[sub_name]] <- sub_dist
  }
}

# saveRDS(dist_list, file.path(wd, "predicted_distribution_by2Char_20251202.rds"))
saveRDS(dist_list, file.path(wd, "predicted_distribution_by2Char_upd2026.rds"))

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

# 
# d4s %>%
#   filter(age70plus == 1, female == 1, site == "galb") %>%
#   group_by(diagnosed_stage) %>%
#   tally() %>%
#   mutate(prop = n / sum(n))


# Use the parameter above for base-case analysis
# use simply recalibration for scenario analysis, by multiply 1.2
# params[, c("tp12", "tp23", "tp34")] <- params[, c("tp12", "tp23", "tp34")]*1.2
# The justification could be the parameters may under-estimate stage 4 cancer,
# therefore underestimate treatment effect

# Re-calibration below may not be needed.

# Re-calibration ----------------------------------------------------------
all_params <- do.call(rbind, params_list)
subgroups <- c(rep("<70; Men", 4), rep("<70; Women", 4), rep(">=70; Men", 4), rep(">=70; Women", 4))
all_params <- cbind(subgroups, all_params)
rownames(all_params) <- NULL

subgroups <- c("<70; Men", "<70; Women", ">=70; Men", ">=70; Women")

recali_list <- list()
for (sub in subgroups){
  
  params <- all_params %>% filter(.data[["subgroups"]]==sub)
  
  char_value1 <- ifelse(grepl(">=70", sub), 1, 0)
  char_value2 <- ifelse(grepl("Women", sub), 1, 0)
  
  sub_data <- d4s %>%
    filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
  
  calibration_results <- iterative_calibration(
    patient_data = sub_data,
    initial_params = params,
    max_iterations = 100,
    convergence_threshold = 0.025,
    adjustment_weight = 0.075,
    n_sim = 200,
    adaptive_weight = TRUE
  )
  
  recali_list[[sub]] <- calibration_results
}
saveRDS(recali_list, file.path(wd, "recali_params_list_v6.rds"))

plot_list_recali <- list()
for (sub in subgroups){
  
  recali_sub <- recali_list[[sub]]
  sub_plot <- plot_comparison(recali_sub$final_comparison, re_cali = T) + ggtitle(sub) + 
    theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "none")
  
  plot_list_recali[[sub]] <- sub_plot
}

# Create one plot with legend to extract it
plot_with_legend <- plot_comparison(recali_list[[sub]]$final_comparison, re_cali = T) + 
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
fig <- grid.arrange(plot_list_recali[[1]], plot_list_recali[[2]], 
                    plot_list_recali[[3]], plot_list_recali[[4]], 
                    legend,
                    ncol = 2, nrow = 3,
                    heights = c(1, 1, 0.2))

ggsave(file.path(output, paste0(opt_tag, "_valiation_by2Char_recali_v6.png")), fig, width = 40, height = 50, units = "cm", dpi = 300)


# v2
recali_params <- do.call(rbind, lapply(recali_list, function(x) x$calibrated_params))
recali_params[, grepl("tp", colnames(recali_params))] <- round(recali_params[, grepl("tp", colnames(recali_params))], 3)
rownames(recali_params) <- NULL
recali_params %>% arrange(site, subgroups) %>% relocate(site)













