# Estimate treatment effect
# Apply to the Cancer progression model
rm(list = ls())
library(tidyverse)
library(survival)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# pop_nam <- "symp_stageImputed"
pop_nam <- "symp_stageImputed_upd2026"
# "symp"
# "noSymp_3rdCons"
# "noSymp_ConsB4act"
# "sympPLUSnoSymp_3rdCons"
# "sympPLUSnoSymp_ConsB4act"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds")))
  
# Check interaction between recommendations and sex and age ---------------
fm <- as.formula("time2diag_surv ~ age70plus + female + nonwhite + imd5_imp2 + ng12_red_flag2 + age70plus:ng12_red_flag2 + female:ng12_red_flag2")

sites <- c("oeso", "stom", "panc", "galb")

mod_list <- list()
for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")

  mod_list[[st]] <- wei
}

oeso <- summary(mod_list$oeso)$table[, c("Value", "Std. Error")]
stom <- summary(mod_list$stom)$table[, c("Value", "Std. Error")]
panc <- summary(mod_list$panc)$table[, c("Value", "Std. Error")]
galb <- summary(mod_list$galb)$table[, c("Value", "Std. Error")]

expo <- round(cbind(oeso, stom, panc, galb),2)
coln <- c("oeso", "", "stom", "", "panc", "", "galb", "")
expo <- rbind(coln,  expo)
# write.csv(expo, file.path(output, "aft_wei_time2diag_checkInt.csv"))
write.csv(expo, file.path(output, "aft_wei_time2diag_checkInt_upd2026.csv"))

# Estimate tx for subgroups by age and sex  -------------------------------
## Extract tx coefficient
extract_tx <- function(model) {
  tx_img <- 1/exp(coef(model)["ng12_red_flag2Imaging"])
  tx_2ww <- 1/exp(coef(model)["ng12_red_flag22 Week Wait"])
  rt <- c(tx_img, tx_2ww)
  names(rt) <- c("imaging", "2ww")
  return(rt)
}

sites <- c("oeso", "stom", "panc", "galb")
fm <- as.formula("time2diag_surv ~ nonwhite + imd5_imp2 + ng12_red_flag2")
tx_mod_list <- list()
tx_list <- list()

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    for (st in sites) {
    
    d4s_sub <- subset(d4s, site == st & age70plus== char_value1 & female == char_value2)

    sub_name <- paste0("age70plus",char_value1, "_", "female", char_value2)
    
    time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
    
    wei <- survreg(fm, d4s_sub, dist = "weibull")
    
    tx_mod_list[[sub_name]][[st]] <- wei
    
    tx_list[[sub_name]][[st]] <- extract_tx(wei)
    }
  }
}

# saveRDS(tx_list, file.path(wd, "tx_list.rds"))
saveRDS(tx_list, file.path(wd, "tx_list_upd2026.rds"))

# Estimate tx by site only ------------------------------------------------
# UPDATE 2026-04: use the intervention effect by cancer site only
# However, keep the by age and sex structure to simplicity and flexibility
# tx will be the same for a certain cancer type across age and sex.
sites <- c("oeso", "stom", "panc", "galb")
fm <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")
tx_mod_list <- list()
tx_list <- list()

for (char_value1 in 0:1) {
  for (char_value2 in 0:1){
    for (st in sites) {
      
      d4s_sub <- subset(d4s, site == st)
      
      sub_name <- paste0("age70plus",char_value1, "_", "female", char_value2)
      
      time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
      
      wei <- survreg(fm, d4s_sub, dist = "weibull")
      
      tx_mod_list[[sub_name]][[st]] <- wei
      
      tx_list[[sub_name]][[st]] <- extract_tx(wei)
    }
  }
}

saveRDS(tx_list, file.path(wd, "tx_list_upd2026_bySiteOnly.rds"))

## Present the model results ----
# Function to extract coefficient and SE from summary table
extract_coef_se <- function(model, coef_name) {
  summ_table <- summary(model)$table
  coef_val <- round(summ_table[coef_name, "Value"], 2)
  se_val <- round(summ_table[coef_name, "Std. Error"], 2)
  return(paste0(coef_val, " (", se_val, ")"))
}

# Get all coefficient names (including intercept and Log(scale))
first_model <- tx_mod_list[[1]][[1]]
coef_names <- rownames(summary(first_model)$table)

# Initialize results lists
results_age0 <- list()
results_age1 <- list()

# Extract for each coefficient
for (coef_name in coef_names) {
  
  # Age < 70 (age70plus == 0)
  row_age0 <- c()
  for (fem in 0:1) {
    for (st in sites) {
      sub_name <- paste0("age70plus0_female", fem)
      model <- tx_mod_list[[sub_name]][[st]]
      row_age0 <- c(row_age0, extract_coef_se(model, coef_name))
    }
  }
  results_age0[[coef_name]] <- row_age0
  
  # Age >= 70 (age70plus == 1)
  row_age1 <- c()
  for (fem in 0:1) {
    for (st in sites) {
      sub_name <- paste0("age70plus1_female", fem)
      model <- tx_mod_list[[sub_name]][[st]]
      row_age1 <- c(row_age1, extract_coef_se(model, coef_name))
    }
  }
  results_age1[[coef_name]] <- row_age1
}

# Create column names
col_names <- c(
  paste0("Male_", sites),
  paste0("Female_", sites)
)

# Convert to data frames
table_age0 <- as.data.frame(do.call(rbind, results_age0))
colnames(table_age0) <- col_names
table_age0 <- cbind(Variable = rownames(table_age0), table_age0)
rownames(table_age0) <- NULL

table_age1 <- as.data.frame(do.call(rbind, results_age1))
colnames(table_age1) <- col_names
table_age1 <- cbind(Variable = rownames(table_age1), table_age1)
rownames(table_age1) <- NULL

# write.csv(table_age0, file.path(output, "aft_wei_age70plus0.csv"), row.names = FALSE)
# write.csv(table_age1, file.path(output, "aft_wei_age70plus1.csv"), row.names = FALSE)

write.csv(table_age0, file.path(output, "aft_wei_age70plus0_upd2026.csv"), row.names = FALSE)
write.csv(table_age1, file.path(output, "aft_wei_age70plus1_upd2026.csv"), row.names = FALSE)
