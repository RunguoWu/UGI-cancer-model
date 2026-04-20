# Use the recalibrated natural history model parameters
# Estimate the impact of tx on cancer stage distributions

rm(list = ls())
library(tidyverse)
library(scales)
library(gridExtra)  # or library(patchwork) if you prefer
library(grid)
library(survival)
library(rstpm2)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))
source(file.path(scr, "natural_history_model", "fn_tp_recalibration.R"))

# Patient data
dt <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) 
d4s <- dt %>% 
  # Only look at those without red flag symptoms according to NG12
  # Because those with red flag symptoms are very likely to be in a different route
  filter(ng12_red_flag == "No red flag") 

# Model parameters
# recali_list <- readRDS(file.path(wd, "recali_params_list_v6.rds"))
# params <- do.call(rbind, lapply(recali_list, function(x) x$calibrated_params))
# rownames(params) <- NULL
# params %>% arrange(site, subgroups) %>% relocate(site)
# For easier access later
# names(recali_list) <- c("age70plus0_female0", "age70plus0_female1", 
#                         "age70plus1_female0", "age70plus1_female1")

params_list <- readRDS(file.path(wd, "params_list_20251202.rds"))

# Load tx
tx_list <- readRDS(file.path(wd, "tx_list.rds"))

# Subgroups
site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
tx_name <- "2ww" # "imaging" #   


# Snr 1: full tx ----------------------------------------------------------
output_list <- list()
ind_output_list <- list()

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    
    sub_name <- paste0(char_name[1],char_value1, "_", char_name[2], char_value2)
    
    sub_data <- d4s %>%
      filter(.data[[char_name[1]]]== char_value1, .data[[char_name[2]]]== char_value2)
    
    # params <- recali_list[[sub_name]][["calibrated_params"]]
    params <- params_list[[sub_name]]
    
    tx <- sapply(tx_list[[sub_name]], "[", tx_name)
    names(tx) <- gsub(paste0("\\.", tx_name, "$"), "", names(tx))
    
    pred_dist_interv <- predict_stage_distribution2(sub_data, params, 
                                                    n_sim = 1000, months=24, 
                                                    tx=tx, tx_prob=1)
    
    comparison_df <- compare_distributions_interv(pred_dist_interv)
    
    output_list[[sub_name]] <- comparison_df
    
    ind_output_list[[sub_name]] <- pred_dist_interv
  }
}
# saveRDS(output_list, file.path(output, paste0("interv_", tx_name, "_by2Char.rds")))
saveRDS(ind_output_list, file.path(output, paste0("ind_interv_", tx_name, "_by2Char.rds")))

# tx_name <- "imaging"
tx_name <- "2ww"
ind_output_list <- readRDS(file.path(output, paste0("ind_interv_", tx_name, "_by2Char.rds")))

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

ggsave(file.path(output, paste0("interv_", tx_name, "_by2Char.png")), fig, width = 28, height = 35, units = "cm", dpi = 300)


# Long-term projection ----------------------------------------------------

# Load intervention result
# tx_name <- "imaging"
# tx_name <- "2ww"

ind_output_list_img <- readRDS(file.path(output, paste0("ind_interv_", "imaging", "_by2Char.rds")))
ind_output_list_2ww <- readRDS(file.path(output, paste0("ind_interv_", "2ww", "_by2Char.rds")))

mod_list <- readRDS(file.path(wd, "stpm2_mod_list.rds" ))

stage_list <- list()
for (st in c("panc", "oeso", "stom", "galb")){
  stage_list[[st]] <- NULL
  for (i in 1:length(ind_output_list_img)) {
    
    e_patid <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["e_patid"]]
    stage_proportions_img <- ind_output_list_img[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
    
    stage_proportions_2ww <- ind_output_list_2ww[[i]][[st]][["patient_level_results"]][["stage_proportions"]]
    
    stage_proportions_df <- data.frame(
      e_patid = e_patid,
      X1_img = stage_proportions_img[, 1],
      X2_img = stage_proportions_img[, 2],
      X3_img = stage_proportions_img[, 3],
      X4_img = stage_proportions_img[, 4],
      X1_2ww = stage_proportions_2ww[, 1],
      X2_2ww = stage_proportions_2ww[, 2],
      X3_2ww = stage_proportions_2ww[, 3],
      X4_2ww = stage_proportions_2ww[, 4],
      stringsAsFactors = FALSE)
    
    stage_list[[st]] <- rbind(stage_list[[st]], stage_proportions_df)
  }
}

stage_df <- do.call(rbind, stage_list)

d4s2 <- d4s %>% 
  select(e_patid, site, female, age10_cent60, age10_new, age70plus, nonwhite,  
         imd5_imp2, stage_imp, death_cancer) %>% 
  left_join(stage_df)

## Descriptive summary----
summary_table <- d4s2 %>%
  group_by(site) %>%
  summarise(
    # Observed proportions at each stage
    obs_stage1 = mean(stage_imp == 1),
    obs_stage2 = mean(stage_imp == 2),
    obs_stage3 = mean(stage_imp == 3),
    obs_stage4 = mean(stage_imp == 4),
    # Mean intervention probabilities
    int_img_stage1 = mean(X1_img),
    int_img_stage2 = mean(X2_img),
    int_img_stage3 = mean(X3_img),
    int_img_stage4 = mean(X4_img),
    
    int_2ww_stage1 = mean(X1_2ww),
    int_2ww_stage2 = mean(X2_2ww),
    int_2ww_stage3 = mean(X3_2ww),
    int_2ww_stage4 = mean(X4_2ww),
    # Total number of patients
    n = n()
  )

summary_long <- summary_table %>%
  pivot_longer(
    cols = -c(site, n),
    names_to = c("type", "stage"),
    names_pattern = "(.*)_stage(.*)",
    values_to = "proportion"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = proportion
  ) %>% 
  mutate(
    dif_img = int_img - obs,
    dif_2ww = int_2ww - obs
  )

## Model prediction ----
surv_years <- 5

pred_list <- list()
for (st in site_name) {
  
  mod <- mod_list[[st]]
  st_data <- d4s2 %>% filter(site == st)
  
  pred_list[[st]] <- list()
  pred_st <- NULL
  
  for (stage in as.character(1:4)) {
    
    newdata_int <- st_data %>% 
      mutate(fu_diag2cens = surv_years,
             stage_imp = stage)
    
    pred_int <- predict(mod, newdata = newdata_int, type = "surv")
    pred_st <- cbind(pred_st, pred_int)
  }
  
  newdata_cur <- st_data %>% mutate(fu_diag2cens = surv_years)
  pred_cur <- predict(mod, newdata = newdata_cur, type = "surv")
  
  surv_st <- data.frame(st_data$e_patid, pred_st, pred_cur)
  colnames(surv_st) <- c("e_patid", paste0("diag_s", 1:4), "surv_cur")
  
  pred_list[[st]] <- surv_st
}

pred_dt <- do.call(rbind, pred_list)

# Merge back to patient data
d4s2 <- d4s2 %>% left_join(pred_dt)

# Calculate survival under intervention
d4s2 <- d4s2 %>% 
  mutate(
    surv_int_img = X1_img * diag_s1 + X2_img * diag_s2 + X3_img * diag_s3 + X4_img * diag_s4,
    surv_int_2ww = X1_2ww * diag_s1 + X2_2ww * diag_s2 + X3_2ww * diag_s3 + X4_2ww * diag_s4
  )

d4s2 %>% group_by(site) %>% 
  summarise(mean(surv_cur), mean(surv_int_img), dif_img = mean(surv_int_img - surv_cur),
            mean(surv_int_2ww), dif_2ww = mean(surv_int_2ww - surv_cur)
            )


d4s2 %>% group_by(site) %>% 
  summarise(mean(surv_cur), mean(surv_int), imp = mean(surv_int - surv_cur))

d4s2 %>% group_by(female) %>% 
  summarise(mean(surv_cur), mean(surv_int), imp = mean(surv_int - surv_cur))

d4s2 %>% group_by(age70plus) %>% 
  summarise(mean(surv_cur), mean(surv_int), imp = mean(surv_int - surv_cur))







