# Fit flexible parametric survival model 
# Predict long-term survival after diagnosis

rm(list = ls())

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "tx_model", "fn_analysis.R"))

library(tidyverse)
library(survival)
library(rstpm2)
library(fastDummies)
library(gridExtra)
library(grid)
library(tidyr)

# pop_nam <- "symp_stageImputed"
pop_nam <- "symp_stageImputed_upd2026"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) %>% 
  mutate(death = if_else(is.na(death_date), 0, 1),
         death_cancer = if_else(is.na(death_cancer), 0, death_cancer),
         death_upGI = if_else(is.na(death_upGI), 0, death_cancer),
         death_oeso = if_else(is.na(death_oeso), 0, death_cancer),
         death_stom = if_else(is.na(death_stom), 0, death_cancer),
         death_panc = if_else(is.na(death_panc), 0, death_cancer),
         death_galb = if_else(is.na(death_galb), 0, death_cancer)
         ) %>% 
  mutate(
    fu_diag2cens = as.numeric(fu_end_date - cancerdate),
    fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 0.5, fu_diag2cens),
    fu_diag2cens = fu_diag2cens/365.25
  ) %>% 
  select(e_patid, site, fu_diag2cens, death_cancer, age10_cent60, nonwhite, female,
           imd5_imp2, stage_imp)


# Model selection ---------------------------------------------------------
mod_list <- list()

## oeso -----
st <- "oeso"
d4s_sub <- subset(d4s, site == st)

fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 + nonwhite + female +
   imd5_imp2 + stage_imp 

fit <- stpm2(fm, d4s_sub, df=1)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=2)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=3)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=4)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=5)
AIC(fit)

fit <- stpm2(fm, d4s_sub, df=5)
plot_stpm2(fit, by_stage = 1, max_time = 5, d4s_sub)

mod_list[[st]] <- fit

## stom -----
st <- "stom"
d4s_sub <- subset(d4s, site == st)

fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 + nonwhite + female +
  imd5_imp2 + stage_imp 

fit <- stpm2(fm, d4s_sub, df=1)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=2)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=3)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=4)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=5)
AIC(fit)

fit <- stpm2(fm, d4s_sub, df=5)
plot_stpm2(fit, by_stage = 1, max_time = 5, d4s_sub)
mod_list[[st]] <- fit

## panc -----
st <- "panc"
d4s_sub <- subset(d4s, site == st)

fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60*stage_imp + nonwhite + female +
  imd5_imp2 + stage_imp 

fit <- stpm2(fm, d4s_sub, df=1)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=2)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=3)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=4)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=5)
AIC(fit)

fit <- stpm2(fm, d4s_sub, df=5)
plot_stpm2(fit, by_stage = 1, max_time = 5, d4s_sub)
mod_list[[st]] <- fit

## galb -----
st <- "galb"
d4s_sub <- subset(d4s, site == st)

fm <- Surv(fu_diag2cens, death_cancer) ~ age10_cent60 + nonwhite + female +
  imd5_imp2 + stage_imp 

fit <- stpm2(fm, d4s_sub, df=1)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=2)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=3)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=4)
AIC(fit)
fit <- stpm2(fm, d4s_sub, df=5)
AIC(fit)

fit <- stpm2(fm, d4s_sub, df=5)
plot_stpm2(fit, by_stage = 1, max_time = 5, d4s_sub)
mod_list[[st]] <- fit

# saveRDS(mod_list, file.path(wd, "stpm2_mod_list.rds" ))
saveRDS(mod_list, file.path(wd, "stpm2_mod_list_upd2026.rds" ))

# After selection ---------------------------------------------------------
mod_list <- readRDS(file.path(wd, "stpm2_mod_list_upd2026.rds" ))

p1 <- plot_stpm2(mod_list[["oeso"]], by_stage=1, max_time=5, 
                        subset(d4s, site=="oeso"), add_legend=FALSE) +
  ggtitle("Oesophagus") + theme(plot.title = element_text(size = 14, face = "bold"))

p2 <- plot_stpm2(mod_list[["stom"]], by_stage=1, max_time=5, 
                 subset(d4s, site=="stom"), add_legend=FALSE) +
  ggtitle("Stomach") + theme(plot.title = element_text(size = 14, face = "bold"))

p3 <- plot_stpm2(mod_list[["panc"]], by_stage=1, max_time=5, 
                 subset(d4s, site=="panc"), add_legend=FALSE) +
  ggtitle("Pancreas") + theme(plot.title = element_text(size = 14, face = "bold"))

p4 <- plot_stpm2(mod_list[["galb"]], by_stage=1, max_time=5, 
                 subset(d4s, site=="galb"), add_legend=FALSE) +
  ggtitle("Gallbladder") + theme(plot.title = element_text(size = 14, face = "bold"))

p_legend <- plot_stpm2(mod_list[["oeso"]], by_stage=1, max_time=5, subset(d4s, site=="oeso"), add_legend=TRUE)

# Extract the legend
get_legend <- function(p){
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p_legend)

# Arrange plots with shared legend
fig <- grid.arrange(
  arrangeGrob(p3, p1, p2, p4, ncol=2),
  legend,
  ncol=2,
  widths=c(8, 1.5)
)

# ggsave(file.path(output, "stpm2_mod.png"), fig, width = 20, height = 18, units = "cm", dpi = 300)
ggsave(file.path(output, "stpm2_mod_upd2026.png"), fig, width = 20, height = 18, units = "cm", dpi = 300)

# Present models ----------------------------------------------------------

# Function to extract coefficients from stpm2 model
extract_coef <- function(model, model_name) {
  coef_summary <- coef(summary(model))
  
  data.frame(
    cancer = model_name,
    parameter = rownames(coef_summary),
    estimate = sprintf("%.2f", coef_summary[, "Estimate"]),
    sd = sprintf("%.2f", coef_summary[, "Std. Error"]),
    stringsAsFactors = FALSE
  )
}

# Extract coefficients from all models
coef_list <- list(
  extract_coef(mod_list$oeso, "Oesophagus"),
  extract_coef(mod_list$stom, "Stomach"),
  extract_coef(mod_list$panc, "Pancreas"),
  extract_coef(mod_list$galb, "Gallbladder")
)

# Combine into one data frame
coef_df <- bind_rows(coef_list)

# Reshape to wide format for easier comparison
coef_wide <- coef_df %>%
  mutate(est_sd = paste0(estimate, " (", sd, ")")) %>%
  select(cancer, parameter, est_sd) %>%
  pivot_wider(names_from = cancer, values_from = est_sd)

# Print the table
print(coef_wide)

# Optional: Save as CSV
# write.csv(coef_wide, file.path(output, "stpm2_mod.csv"), row.names = FALSE)
write.csv(coef_wide, file.path(output, "stpm2_mod_upd2026.csv"), row.names = FALSE)




