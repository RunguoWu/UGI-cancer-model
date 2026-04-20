rm(list = ls())
library(tidyverse)
library(survival)
library(gtsummary)
library(skimr)
library(dplyr)
library(rstpm2)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "tx_model", "fn_analysis.R"))

source(file.path(scr, "natural_history_model", "fn_tp_optimise_4stage.R"))
source(file.path(scr, "natural_history_model", "fn_parameter_search.R"))
source(file.path(scr, "natural_history_model", "fn_tp_validation.R"))

# pop_nam <- "symp_stageImputed_1yrIndex"
pop_nam <- "symp_stageImputed_1yrIndex_upd2026"
d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) 

# Table 1 -----------------------------------------------------------------

# keep missing ethnicity for summary only
d4s <- d4s %>% 
  mutate(
    ethn2_mis = if_else(is.na(ethnicity), NA, ethn2)
  )

rt <- d4s %>% select(age_index, gender, ethn2_mis, imd5_origin, site, stage,
                     ng12_red_flag2, action_type22, 
                     time2diag, time2act, time_act2diag) %>%
  tbl_summary(missing = "ifany",
              statistic = list(all_continuous() ~ "{mean} ({sd})")
  ) %>%
  add_n() %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = file.path(output, paste0("table1_", pop_nam, ".docx")))


# Table 2 -----------------------------------------------------------------

## Overall One-step model -------------------------------------------------

fm <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")

sites <- c("oeso", "stom", "panc", "galb")

mod_list <- list()
for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")
  exp <- survreg(fm, d4s_sub, dist = "exponential")
  lgn <- survreg(fm, d4s_sub, dist = "lognormal")
  lgl <- survreg(fm, d4s_sub, dist = "loglogistic")
  
  mod_list[[st]][["wei"]] <- wei
  mod_list[[st]][["exp"]] <- exp
  mod_list[[st]][["lgn"]] <- lgn
  mod_list[[st]][["lgl"]] <- lgl
}

# interactions between red flags and other characteristics have been checked
# no very strong interactions
# keep the current one
######################################################################
AIC_list <- data.frame(
  "oeso" = sapply(mod_list[["oeso"]], AIC),
  "stom" = sapply(mod_list[["stom"]], AIC),
  "panc" = sapply(mod_list[["panc"]], AIC),
  "galb" = sapply(mod_list[["galb"]], AIC)
)
######################################################################
# AIC suggest Weibull for all
oeso <- summary(mod_list$oeso$wei)$table[, c("Value", "Std. Error")]
stom <- summary(mod_list$stom$wei)$table[, c("Value", "Std. Error")]
panc <- summary(mod_list$panc$wei)$table[, c("Value", "Std. Error")]
galb <- summary(mod_list$galb$wei)$table[, c("Value", "Std. Error")]

oeso1 <- paste0(round(oeso[, "Value"],2), " (", round(oeso[, "Std. Error"],2), ")")
stom1 <- paste0(round(stom[, "Value"],2), " (", round(stom[, "Std. Error"],2), ")")
panc1 <- paste0(round(panc[, "Value"],2), " (", round(panc[, "Std. Error"],2), ")")
galb1 <- paste0(round(galb[, "Value"],2), " (", round(galb[, "Std. Error"],2), ")")

expo <- cbind(panc1, oeso1, stom1, galb1)
coln <- c("Pancreatic Weibull", "Oesophageal Weibull", "Gastric Weibull", 
          "Gallbladder Weibull")
expo <- rbind(coln,  expo)
colnames(expo) <- NULL
rownames(expo) <- c("Variables", "Intercept", "60-69 years", "70-79 years",
                    "80 years plus", "Women", "Non-White", "IMD Q2", "IMD Q3",
                    "IMD Q4", "IMD Q5", "Rcmd Imaging", "Rcmd cancer referral",
                    "Log(scale)"
)
dt <- as.data.frame(expo[-1, ])
colnames(dt) <- expo[1, ]

# write.csv(dt, file.path(output, "aft_1step_1yrIndex.csv"))
write.csv(dt, file.path(output, "aft_1step_1yrIndex_upd2026.csv"))

## Two-step model ----------------------------------------------------------
### time to action ----
fm <- as.formula("time2act_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")

sites <- c("panc", "oeso", "stom", "galb")

mod_list <- list()
for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  d4s_sub$time2act <- ifelse(d4s_sub$time2act==0, 1, d4s_sub$time2act)
  
  time2act_surv <- Surv(time = d4s_sub$time2act, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")
  exp <- survreg(fm, d4s_sub, dist = "exponential")
  lgn <- survreg(fm, d4s_sub, dist = "lognormal")
  lgl <- survreg(fm, d4s_sub, dist = "loglogistic")
  
  mod_list[[st]][["wei"]] <- wei
  mod_list[[st]][["exp"]] <- exp
  mod_list[[st]][["lgn"]] <- lgn
  mod_list[[st]][["lgl"]] <- lgl
}

######################################################################
AIC_list <- data.frame(
  "oeso" = sapply(mod_list[["oeso"]], AIC),
  "stom" = sapply(mod_list[["stom"]], AIC),
  "panc" = sapply(mod_list[["panc"]], AIC),
  "galb" = sapply(mod_list[["galb"]], AIC)
)
######################################################################
oeso <- summary(mod_list$oeso$wei)$table[, c("Value", "Std. Error")]
stom <- summary(mod_list$stom$wei)$table[, c("Value", "Std. Error")]
panc <- summary(mod_list$panc$wei)$table[, c("Value", "Std. Error")]
galb <- summary(mod_list$galb$wei)$table[, c("Value", "Std. Error")]

oeso1 <- paste0(round(oeso[, "Value"],2), " (", round(oeso[, "Std. Error"],2), ")")
stom1 <- paste0(round(stom[, "Value"],2), " (", round(stom[, "Std. Error"],2), ")")
panc1 <- paste0(round(panc[, "Value"],2), " (", round(panc[, "Std. Error"],2), ")")
galb1 <- paste0(round(galb[, "Value"],2), " (", round(galb[, "Std. Error"],2), ")")

expo_p1 <- cbind(panc1, oeso1, stom1, galb1)
coln <- c("Pancreatic P1", "Oesophageal P1", "Gastric P1", 
          "Gallbladder P1")
expo_p1 <- rbind(coln,  expo_p1)
colnames(expo_p1) <- NULL
rownames(expo_p1) <- c("Variables", "Intercept", "60-69 years", "70-79 years",
                       "80 years plus", "Women", "Non-White", "IMD Q2", "IMD Q3",
                       "IMD Q4", "IMD Q5", "Rcmd Imaging", 
                       "Rcmd cancer referral",
                       "Log(scale)"
)

### time from action to diagnosis ----
# check effect in regression coefficients by stage and sites
fm <- as.formula("time_act2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + action_type22")
sites <- c("panc", "oeso", "stom", "galb")

mod_list <- list()
for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  d4s_sub$time_act2diag <- ifelse(d4s_sub$time_act2diag==0, 1, d4s_sub$time_act2diag)
  
  time_act2diag_surv <- Surv(time = d4s_sub$time_act2diag, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")
  exp <- survreg(fm, d4s_sub, dist = "exponential")
  lgn <- survreg(fm, d4s_sub, dist = "lognormal")
  lgl <- survreg(fm, d4s_sub, dist = "loglogistic")
  
  mod_list[[st]][["wei"]] <- wei
  mod_list[[st]][["exp"]] <- exp
  mod_list[[st]][["lgn"]] <- lgn
  mod_list[[st]][["lgl"]] <- lgl
}

######################################################################
AIC_list <- data.frame(
  "oeso" = sapply(mod_list[["oeso"]], AIC),
  "stom" = sapply(mod_list[["stom"]], AIC),
  "panc" = sapply(mod_list[["panc"]], AIC),
  "galb" = sapply(mod_list[["galb"]], AIC)
)
# Use Weibull for all, checked by plotting
# lgn gives some unexplanable results for oeso and stom.

######################################################################
oeso <- summary(mod_list$oeso$wei)$table[, c("Value", "Std. Error")]
stom <- summary(mod_list$stom$wei)$table[, c("Value", "Std. Error")]
panc <- summary(mod_list$panc$wei)$table[, c("Value", "Std. Error")]
galb <- summary(mod_list$galb$wei)$table[, c("Value", "Std. Error")]

oeso1 <- paste0(round(oeso[, "Value"],2), " (", round(oeso[, "Std. Error"],2), ")")
stom1 <- paste0(round(stom[, "Value"],2), " (", round(stom[, "Std. Error"],2), ")")
panc1 <- paste0(round(panc[, "Value"],2), " (", round(panc[, "Std. Error"],2), ")")
galb1 <- paste0(round(galb[, "Value"],2), " (", round(galb[, "Std. Error"],2), ")")

expo_p2 <- cbind(panc1, oeso1, stom1, galb1)
coln <- c("Pancreatic P2", "Oesophageal P2", "Gastric P2", "Gallbladder P2")
expo_p2 <- rbind(coln,  expo_p2)
colnames(expo_p2) <- NULL
rownames(expo_p2) <- c("Variables", "Intercept", "60-69 years", "70-79 years",
                       "80 years plus", "Women", "Non-White", "IMD Q2", "IMD Q3",
                       "IMD Q4", "IMD Q5", "Imaging", 
                       "Cancer referral",
                       "Log(scale)"
)

dt_p1 <- as.data.frame(expo_p1[-1, ])
colnames(dt_p1) <- expo_p1[1, ]

dt_p2 <- as.data.frame(expo_p2[-1, ])
colnames(dt_p2) <- expo_p2[1, ]

dt_combined <- merge(dt_p1, dt_p2, by = "row.names", all = TRUE)
rownames(dt_combined) <- dt_combined$Row.names
dt_combined$Row.names <- NULL

new_order <- c("Intercept", "60-69 years", "70-79 years",
               "80 years plus", "Women", "Non-White", "IMD Q2", "IMD Q3",
               "IMD Q4", "IMD Q5", "Rcmd Imaging", 
               "Rcmd cancer referral", "Imaging", 
               "Cancer referral", "Log(scale)")
dt_combined <- dt_combined[new_order, c("Pancreatic P1", "Pancreatic P2",
                                        "Oesophageal P1", "Oesophageal P2",
                                        "Gastric P1", "Gastric P2",
                                        "Gallbladder P1", "Gallbladder P2"
)]

# write.csv(dt_combined, file.path(output, "aft_2step_1yrIndex.csv"))
write.csv(dt_combined, file.path(output, "aft_2step_1yrIndex_upd2026.csv"))

dt_3 <- merge(dt, dt_combined, by = "row.names", all = TRUE)
rownames(dt_3) <- dt_3$Row.names
dt_3$Row.names <- NULL

dt_3 <- dt_3[new_order, ]
# write.csv(dt_3, file.path(output, "aft_1&2step_1yrIndex.csv"))
write.csv(dt_3, file.path(output, "aft_1&2step_1yrIndex_upd2026.csv"))


# Table 3 -----------------------------------------------------------------
## simulate the full tx first ----
 
# d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_1yrIndex.rds")) %>%
d4s <- readRDS(file.path(wd, "study_pop_symp_stageImputed_1yrIndex_upd2026.rds")) %>%
  filter(ng12_red_flag == "No red flag") 

# params_list <- readRDS(file.path(wd, "params_list_1yrIndex.rds"))
params_list <- readRDS(file.path(wd, "params_list_1yrIndex_upd2026.rds"))

# Load tx
# tx_list <- readRDS(file.path(wd, "tx_list_1yrIndex.rds"))
# tx_list <- readRDS(file.path(wd, "tx_list_1yrIndex_upd2026.rds"))
tx_list <- readRDS(file.path(wd, "tx_list_1yrIndex_upd2026_bySiteOnly.rds"))

# Subgroups
site_name <- c("galb", "oeso", "panc", "stom")
char_name <- c("age70plus", "female")
tx_name <- "imaging" # "2ww" #    

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
                                                    n_sim = 1000, months=12, 
                                                    tx=tx, tx_prob=1)
    
    comparison_df <- compare_distributions_interv(pred_dist_interv)
    
    output_list[[sub_name]] <- comparison_df
    
    ind_output_list[[sub_name]] <- pred_dist_interv
  }
}
# saveRDS(output_list, file.path(output, paste0("interv_", tx_name, "_by2Char_1yrIndex.rds")))
# saveRDS(ind_output_list, file.path(output, paste0("ind_interv_", tx_name, "_by2Char_1yrIndex.rds")))
# saveRDS(output_list, file.path(output, paste0("interv_", tx_name, "_by2Char_1yrIndex_upd2026.rds")))
# saveRDS(ind_output_list, file.path(output, paste0("ind_interv_", tx_name, "_by2Char_1yrIndex_upd2026.rds")))
saveRDS(output_list, file.path(output, paste0("interv_", tx_name, "_by2Char_1yrIndex_upd2026_bySiteOnly.rds")))
saveRDS(ind_output_list, file.path(output, paste0("ind_interv_", tx_name, "_by2Char_1yrIndex_upd2026_bySiteOnly.rds")))

## Scenario analyses
# ind_output_list_img <- readRDS(file.path(output, paste0("ind_interv_", "imaging", "_by2Char_1yrIndex.rds")))
# ind_output_list_2ww <- readRDS(file.path(output, paste0("ind_interv_", "2ww", "_by2Char_1yrIndex.rds")))
# ind_output_list_img <- readRDS(file.path(output, paste0("ind_interv_", "imaging", "_by2Char_1yrIndex_upd2026.rds")))
# ind_output_list_2ww <- readRDS(file.path(output, paste0("ind_interv_", "2ww", "_by2Char_1yrIndex_upd2026.rds")))
ind_output_list_img <- readRDS(file.path(output, paste0("ind_interv_", "imaging", "_by2Char_1yrIndex_upd2026_bySiteOnly.rds")))
ind_output_list_2ww <- readRDS(file.path(output, paste0("ind_interv_", "2ww", "_by2Char_1yrIndex_upd2026_bySiteOnly.rds")))

# mod_list <- readRDS(file.path(wd, "stpm2_mod_list.rds" ))
mod_list <- readRDS(file.path(wd, "stpm2_mod_list_upd2026.rds" ))

site_name <- c("panc", "oeso", "stom", "galb")

## predict stage shift ----
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

# Merge back to patient data for analysis
d4s2 <- d4s %>% 
  select(e_patid, site, female, age10_cent60, age10_new, age70plus, nonwhite,  
         imd5_imp2, stage_imp, death_cancer) %>% 
  left_join(stage_df)

## predict survival ----
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

d4s2 <- d4s2 %>% 
  left_join(pred_dt)

## Scenarios ----
s_img_all <- partial_intervention_summary(d4s2, 
                                          pct_img = 1, 
                                          pct_2ww = 0,
                                          seed = 123)

s_2ww_all <- partial_intervention_summary(d4s2, 
                                          pct_img = 0, 
                                          pct_2ww = 1,
                                          seed = 123)

s_img_haf <- partial_intervention_summary(d4s2, 
                                          pct_img = 0.5, 
                                          pct_2ww = 0,
                                          seed = 123)

s_2ww_haf <- partial_intervention_summary(d4s2, 
                                          pct_img = 0, 
                                          pct_2ww = 0.5,
                                          seed = 123)

s_img_2ww_qua <- partial_intervention_summary(d4s2, 
                                              pct_img = 0.25, 
                                              pct_2ww = 0.25,
                                              seed = 123)

# rt <- cbind(s_img_all[, c("site", "stage", "n_total", "obs_stage", "obs_surv", 
#                           "pred_stage")], 
#             s_2ww_all[, c("pred_stage")],
#             s_img_haf[, c("pred_stage")],
#             s_2ww_haf[, c("pred_stage")],
#             s_img_2ww_qua[, c("pred_stage")],
#             s_img_all[, c("pred_surv")],
#             s_2ww_all[, c("pred_surv")],
#             s_img_haf[, c("pred_surv")],
#             s_2ww_haf[, c("pred_surv")],
#             s_img_2ww_qua[, c("pred_surv")]
# )
# 
# rt[, 4:15] <- round(rt[, 4:15], 3)
# colnames(rt) <- c("Site", "Stage", "N", "Stage distribution", "Five-year survival", 
#                   "pred_stage_img_all",  
#                   "pred_stage_2ww_all", 
#                   "pred_stage_img_haf", 
#                   "pred_stage_2ww_haf", 
#                   "pred_stage_img_2ww_qua", 
#                   "pred_surv_img_all", 
#                   "pred_surv_2ww_all",
#                   "pred_surv_img_haf",
#                   "pred_surv_2ww_haf",
#                   "pred_surv_img_2ww_qua"
# )

# re-organise the summary table
rt <- cbind(s_img_all[, c("site", "stage", "n_total", "obs_stage", "obs_surv", 
                          "pred_stage")], 
            s_img_all[, c("pred_surv")],
            s_2ww_all[, c("pred_stage")],
            s_2ww_all[, c("pred_surv")],
            s_img_haf[, c("pred_stage")],
            s_img_haf[, c("pred_surv")],
            s_2ww_haf[, c("pred_stage")],
            s_2ww_haf[, c("pred_surv")],
            s_img_2ww_qua[, c("pred_stage")],
            s_img_2ww_qua[, c("pred_surv")]
)

rt[, 4:15] <- round(rt[, 4:15], 3) * 100
colnames(rt) <- c("Site", "Stage", "N", "Stage distribution", "Five-year survival", 
                  "pred_stage_img_all",  
                  "pred_surv_img_all", 
                  "pred_stage_2ww_all",
                  "pred_surv_2ww_all", 
                  "pred_stage_img_haf", 
                  "pred_surv_img_haf",
                  "pred_stage_2ww_haf",
                  "pred_surv_2ww_haf",
                  "pred_stage_img_2ww_qua", 
                  "pred_surv_img_2ww_qua"
)


# write.csv(rt, file.path(output, "pred_stage_surv_1yrIndex.csv"))
# write.csv(rt, file.path(output, "pred_stage_surv_1yrIndex_upd2026.csv"))
write.csv(rt, file.path(output, "pred_stage_surv_1yrIndex_upd2026_bySiteOnly.csv"))

