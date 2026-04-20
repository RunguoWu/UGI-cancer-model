# Descriptive analysis of the cancer cohort
# Also brief exploratory analysis 

rm(list = ls())
library(tidyverse)
library(survival)
library(gtsummary)
library(skimr)
library(dplyr)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Descriptive -------------------------------------------------------------
pop_nam <- "symp_stageImputed"
# "symp"
# "noSymp_3rdCons"
# "noSymp_ConsB4act"
# "sympPLUSnoSymp_3rdCons"
# "sympPLUSnoSymp_ConsB4act"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) 

# skim(d4s)

# Tab ---------------------------------------------------------------------

table(d4s$ng12_red_flag2, d4s$action_type21, useNA = "ifany")


# One step model ----------------------------------------------------------
fm <- as.formula("time2diag ~ age70plus + male + ethn2 + imd5_imp2 + ng12_red_flag2")

fit_oeso <- lm(fm, subset(d4s, site == "oeso"))
fit_stom <- lm(fm, subset(d4s, site == "stom"))
fit_panc <- lm(fm, subset(d4s, site == "panc"))
fit_galb <- lm(fm, subset(d4s, site == "galb"))

# m1 <- summary(fit_s1)$coefficients[, c("Estimate", "Std. Error")]
# m2 <- summary(fit_s2)$coefficients[, c("Estimate", "Std. Error")]
# m3 <- summary(fit_s3)$coefficients[, c("Estimate", "Std. Error")]
# m4 <- summary(fit_s4)$coefficients[, c("Estimate", "Std. Error")]
m5 <- summary(fit_oeso)$coefficients[, c("Estimate", "Std. Error")]
m6 <- summary(fit_stom)$coefficients[, c("Estimate", "Std. Error")]
m7 <- summary(fit_panc)$coefficients[, c("Estimate", "Std. Error")]
m8 <- summary(fit_galb)$coefficients[, c("Estimate", "Std. Error")]

expo <- round(cbind(m5, m6, m7, m8),1)
coln <- c("oeso", "", "stom", "", "panc", "", "galb", "")
expo <- rbind(coln,  expo)
write.csv(expo, file.path(output, "ols_time2diag.csv"))


fit <- glm(time2diag ~ action_type4, data=d4s, family = Gamma(link = "log"))
summary(fit)

# try survival AFT
library(survival)
dtt <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds")) %>% 
  filter(ng12_red_flag == "No red flag")
d4s_sub <- dtt # subset(dtt, site == "stom")
time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
fm <- as.formula("time2diag_surv ~ age70plus + female + nonwhite")
aft <- survreg(fm, d4s_sub, dist = "weibull")
summary(aft)

# Two-step model ----------------------------------------------------------

# Boby suggest model time to action and time from action to diagnosis separately,
# Because the action type was not known when symptoms presented so should not be a baseline variable

## time to action ----
fm <- as.formula("time2act ~ age70plus + male + ethn2 + imd5_imp2 + ng12_red_flag2")

# fit_s1 <- lm(fm, subset(d4s, stage == "1"))
# fit_s2 <- lm(fm, subset(d4s, stage == "2"))
# fit_s3 <- lm(fm, subset(d4s, stage == "3"))
# fit_s4 <- lm(fm, subset(d4s, stage == "4"))

fit_oeso <- lm(fm, subset(d4s, site == "oeso"))
fit_stom <- lm(fm, subset(d4s, site == "stom"))
fit_panc <- lm(fm, subset(d4s, site == "panc"))
fit_galb <- lm(fm, subset(d4s, site == "galb"))

# m1 <- summary(fit_s1)$coefficients[, c("Estimate", "Std. Error")]
# m2 <- summary(fit_s2)$coefficients[, c("Estimate", "Std. Error")]
# m3 <- summary(fit_s3)$coefficients[, c("Estimate", "Std. Error")]
# m4 <- summary(fit_s4)$coefficients[, c("Estimate", "Std. Error")]
m5 <- summary(fit_oeso)$coefficients[, c("Estimate", "Std. Error")]
m6 <- summary(fit_stom)$coefficients[, c("Estimate", "Std. Error")]
m7 <- summary(fit_panc)$coefficients[, c("Estimate", "Std. Error")]
m8 <- summary(fit_galb)$coefficients[, c("Estimate", "Std. Error")]

expo <- round(cbind(m5, m6, m7, m8),1)
coln <- c("oeso", "", "stom", "", "panc", "", "galb", "")
expo <- rbind(coln,  expo)
write.csv(expo, file.path(output, "ols_time2action.csv"))

## time from action to diagnosis ----
# check effect in regression coefficients by stage and sites
fm <- as.formula("time_act2diag ~ age70plus + male + ethn2 + imd5_imp2 + action_type21")

# fit_s1 <- lm(fm, subset(d4s, stage == "1"))
# fit_s2 <- lm(fm, subset(d4s, stage == "2"))
# fit_s3 <- lm(fm, subset(d4s, stage == "3"))
# fit_s4 <- lm(fm, subset(d4s, stage == "4"))

fit_oeso <- lm(fm, subset(d4s, site == "oeso"))
fit_stom <- lm(fm, subset(d4s, site == "stom"))
fit_panc <- lm(fm, subset(d4s, site == "panc"))
fit_galb <- lm(fm, subset(d4s, site == "galb"))

# m1 <- summary(fit_s1)$coefficients[, c("Estimate", "Std. Error")]
# m2 <- summary(fit_s2)$coefficients[, c("Estimate", "Std. Error")]
# m3 <- summary(fit_s3)$coefficients[, c("Estimate", "Std. Error")]
# m4 <- summary(fit_s4)$coefficients[, c("Estimate", "Std. Error")]
m5 <- summary(fit_oeso)$coefficients[, c("Estimate", "Std. Error")]
m6 <- summary(fit_stom)$coefficients[, c("Estimate", "Std. Error")]
m7 <- summary(fit_panc)$coefficients[, c("Estimate", "Std. Error")]
m8 <- summary(fit_galb)$coefficients[, c("Estimate", "Std. Error")]

expo <- round(cbind(m5, m6, m7, m8),1)
coln <- c("oeso", "", "stom", "", "panc", "", "galb", "")
expo <- rbind(coln,  expo)

write.csv(expo, file.path(output, "ols_action2diagnosis.csv"))



