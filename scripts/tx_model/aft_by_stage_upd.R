# Estimate treatment effect
# Apply to the Cancer progression model
rm(list = ls())
library(tidyverse)
library(survival)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# pop_nam <- "symp_stageImputed"
pop_nam <- "symp_stageImputed_upd2026"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds")))

d4s$stage_int <- as.integer(d4s$stage_imp)

fm <- as.formula("time2diag_surv ~ age10_new + female + nonwhite + imd5_imp2 + stage_int*ng12_red_flag2")

sites <- c("oeso", "stom", "panc", "galb")

mod_list <- list()
for (st in sites) {
  
  d4s_sub <- subset(d4s, site == st)
  
  time2diag_surv <- Surv(time = d4s_sub$time2diag, event = rep(1, nrow(d4s_sub)))
  
  wei <- survreg(fm, d4s_sub, dist = "weibull")
  mod_list[[st]][["wei"]] <- wei
}

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
                    "IMD Q4", "IMD Q5", "stage as count", "Rcmd Imaging", "Rcmd cancer referral", "stage*Rcmd Imaging", "stage*Rcmd cancer referral", "Log(scale)"
)
dt <- as.data.frame(expo[-1, ])
colnames(dt) <- expo[1, ]

write.csv(dt, file.path(output, paste0("aft_1step", "_stage_upd2026_new.csv")))























d4s %>% group_by(action_type22) %>% summarise(mean(time2diag), mean(time_act2diag))

d4s %>% group_by(ng12_red_flag2, stage_imp) %>% tally()

round(prop.table(table(d4s$ng12_red_flag2, d4s$stage_imp),1),2)

d4s %>% filter(stage_imp=="1") %>% group_by(ng12_red_flag2) %>% summarise(mean(time2diag))
