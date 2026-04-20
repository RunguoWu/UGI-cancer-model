# Estimate treatment effect
# Apply to the Cancer progression model
rm(list = ls())
library(tidyverse)
library(survival)
source("/data/WIPH-CanDetect/HealthEco/route.R")

pop_nam <- "symp_stageImputed_upd2026"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) %>% 
  filter(age_index >=40)

# Estimate tx for subgroups by age and sex  -------------------------------

## Extract tx coefficient
extract_tx <- function(model) {
  tx_img <- 1/exp(coef(model)["ng12_red_flag2Imaging"])
  tx_2ww <- 1/exp(coef(model)["ng12_red_flag22 Week Wait"])
  rt <- c(tx_img, tx_2ww)
  names(rt) <- c("imaging", "2ww")
  return(rt)
}

sites <- c("panc", "oeso", "stom", "galb")
fm <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + ng12_red_flag2")
tx_mod_list <- list()
tx_list <- list()

for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")
  
  tx_mod_list[[st]] <- wei
  
  tx_list[[st]] <- extract_tx(wei)
}


## Present the model results ----
# Function to extract coefficient and SE from summary table
extract_coef_se <- function(model, coef_name) {
  summ_table <- summary(model)$table
  coef_val <- round(summ_table[coef_name, "Value"], 2)
  se_val <- round(summ_table[coef_name, "Std. Error"], 2)
  return(paste0(coef_val, " (", se_val, ")"))
}

# Get all coefficient names (including intercept and Log(scale))
first_model <- tx_mod_list[[1]]
coef_names <- rownames(summary(first_model)$table)

# Initialize results lists
results_age0 <- list()

# Extract for each coefficient
for (coef_name in coef_names) {
  
  # Age < 70 (age70plus == 0)
  row_age0 <- c()
  for (st in sites) {
    model <- tx_mod_list[[st]]
    row_age0 <- c(row_age0, extract_coef_se(model, coef_name))
  }
  
  results_age0[[coef_name]] <- row_age0
}


# Convert to data frames
table_age0 <- as.data.frame(do.call(rbind, results_age0))
table_age0 <- cbind(Variable = rownames(table_age0), table_age0)
rownames(table_age0) <- NULL
colnames(table_age0) <- c("var", sites)
write.csv(table_age0, file.path(output, "aft_wei_age70plus0_upd2026_sa4070.csv"), row.names = FALSE)

