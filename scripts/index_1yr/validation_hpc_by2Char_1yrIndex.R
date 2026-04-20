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

# d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_1yrIndex.rds")) %>% 
d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_1yrIndex_upd2026.rds")) %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") %>%
  select(e_patid, month, diagnosed_stage, site, age70plus, female, nonwhite)

site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
# opt_tag <- "opt_1yrIndex_RS_500K"
opt_tag <- "opt_1yrIndex_upd2026_RS_500K"

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

ggsave(file.path(output, paste0(opt_tag, "convergence2Char_4x4.png")), final_fig, width = 12, height = 10)


# prepare for validation --------------------------------------------------

## Old method ----
# Use the parameters with the maximum likelihood


## New method ----
# Use the median of the distribution of the parameters with likelihood above the threshold
params_mat <- matrix(NA, nrow = 16, ncol=7)
df_list <- list()
rownames(params_mat) <- names(dt_list)
for (sub_name in names(dt_list)) {
  
  dt_list_sub <- dt_list[[sub_name]]
  
  params <- dt_list_sub$final_params
  
  records <- dt_list_sub$params_record
  llh_records <- sapply(records, "[[", "value_cur")
  par_records <- sapply(records, "[[", "par_cur")
  threshold <- quantile(llh_records, probs = 0.01) # .002 corresponding to 1000, with 500K iterations
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
saveRDS(df_list, file.path(wd, "params_list_upd2026_withCI_1yrIndex.rds"))

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

write_xlsx(df_expo, file.path(output, "ugi_cpdi_param_CI_upd2026_1yrIndex.xlsx"))


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

# saveRDS(params_list, file.path(wd, "params_list_1yrIndex.rds"))
saveRDS(params_list, file.path(wd, "params_list_1yrIndex_upd2026.rds"))

# Present parameters
all_params <- do.call(rbind, params_list)
all_params[, grepl("tp", colnames(all_params))] <- round(all_params[, grepl("tp", colnames(all_params))], 3)
subgroups <- c(rep("<70; Men", 4), rep("<70; Women", 4), rep(">=70; Men", 4), rep(">=70; Women", 4))
all_params <- cbind(subgroups, all_params)
rownames(all_params) <- NULL
x <- all_params %>% arrange(site, subgroups) %>% relocate(site)
# write.csv(x, file.path(output, "ugi-cpdi model parameters_1yrIndex.csv"))
write.csv(x, file.path(output, "ugi-cpdi model parameters_1yrIndex_upd2026.csv"))

# Simulate for validation -------------------------------------------------
# Simulate initial stage distribution by site and subgroups
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
              # file.path(output, paste0("initial_dist_1yrIndex", sub_name, ".csv")))
              file.path(output, paste0("initial_dist_1yrIndex_upd2026", sub_name, ".csv")))
  }
}


dist_list <- list()
for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    params <- params_list[[sub_name]]

    sub_dist <- predict_stage_distribution(sub_data,
                                           params, n_sim = 500, months=12,
                                           plus_param = FALSE)
    
    dist_list[[sub_name]] <- sub_dist
  }
}

# saveRDS(dist_list, file.path(wd, "predicted_distribution_by2Char_1yrIndex.rds"))
# dist_list <- readRDS(file.path(wd, "predicted_distribution_by2Char_1yrIndex.rds"))

saveRDS(dist_list, file.path(wd, "predicted_distribution_by2Char_1yrIndex_upd2026.rds"))
dist_list <- readRDS(file.path(wd, "predicted_distribution_by2Char_1yrIndex_upd2026.rds"))

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

# ggsave(file.path(output, paste0(opt_tag, "_valiation_by2Char_1yrIndex.png")), 
#        fig, width = 30, height = 30, units = "cm", dpi = 300)
ggsave(file.path(output, paste0(opt_tag, "_valiation_by2Char_1yrIndex_upd2026.png")), 
       fig, width = 30, height = 30, units = "cm", dpi = 300)

# 
# d4s %>%
#   filter(age70plus == 1, female == 1, site == "galb") %>%
#   group_by(diagnosed_stage) %>%
#   tally() %>%
#   mutate(prop = n / sum(n))








