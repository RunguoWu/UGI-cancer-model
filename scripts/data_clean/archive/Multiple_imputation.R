# Multiply impute missing stage 

rm(list = ls())
library(mice)
library(miceadds)
library(tidyverse)
library(parallel)
library(matrixStats)
library(skimr)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Load data to impute
# dt <- readRDS(file.path(wd, "study_pop_symp.rds"))
dt <- readRDS(file.path(wd, "study_pop_symp_upd2026.rds"))
# Load data for multisite and smoking information
ugi_icd_stage <- readRDS(file.path(wd, "create_study_pop_interm", "ugi_icd_stage.rds"))

dt <-dt %>% left_join(ugi_icd_stage[,c("e_patid", "smokingstatus", "multisite")]) %>% 
  mutate(multisite = if_else(is.na(multisite), 0, multisite))

# Prepare data for impute
dt_imp <- dt %>% 
  mutate(time_diag2end = as.numeric(fu_end_date - cancerdate)) %>% 
  mutate(smoke = if_else(smokingstatus %in% c("current or ex-smoker", "ex-smoker", 
                                              "current smoker", "nicotine or tobacco use",
                                               ""), 1,0)) %>% 
  select(e_patid, gender, age_index, ethn2, imd5_imp2, smoke,
         stage, site, multisite, total_symptoms_28day, risk_level,
         action_type21, time2act, time_act2diag, 
         death, time_diag2end, death_cancer) %>% 
  mutate(
    death_cancer = if_else(is.na(death_cancer), 0, death_cancer),
    stage = if_else(stage == "", NA, stage),
    death_cancer = factor(death_cancer, ordered = F),
    gender = factor(gender, ordered = F),
    ethn2 = factor(ethn2, ordered = F),
    imd5_imp2 = factor(imd5_imp2, levels = c("1", "2", "3", "4", "5"), ordered = T),
    smoke = factor(smoke, ordered = F),
    stage = factor(stage, levels = c("1", "2", "3", "4"), ordered = T),
    site = factor(site, ordered = F),
    multisite = factor(multisite, ordered = F),
    total_symptoms_28day = factor(total_symptoms_28day, ordered = T),
    risk_level = factor(risk_level, ordered = F),
    action_type21 = factor(action_type21, ordered = F),
    death = factor(death, ordered = F)
  )
  
# pre-define the matrix ----
imp <- mice(dt_imp, maxit=0)
predM = imp$predictorMatrix
meth = imp$method
  
# ID and ethn4 won't predict other var.
predM[, c("e_patid")] <- 0 # do not used as predictor
predM[c("e_patid"), ] <- 0 # don't impute 

  
# start impute ----
N <- 20
rand_seed=12345

imp <- futuremice(data = dt_imp, m=N, n.core = 4, parallelseed = rand_seed,
                  predictorMatrix = predM, method = meth)

imp_list <- miceadds::mids2datlist(imp)

saveRDS(imp, file.path(wd, "imp.rds"))
saveRDS(imp_list, file.path(wd, "imp_list.rds"))

# use the mode of imputed values across the 20 imputations as the imputed value ---- 
# Step 1: Combine all imputed datasets with imputation number
# imp_list <- readRDS(file.path(wd, "imp_list.rds"))

combined_imputations <- bind_rows(
  lapply(1:length(imp_list), function(i) {
    imp_list[[i]] %>% 
      mutate(imp_num = i) %>%
      select(e_patid, stage, imp_num)  # Only keep relevant columns
  })
)

# Step 2: Identify patients with originally missing stage values
missing_patients <- dt_imp %>%
  filter(is.na(stage)) %>%
  pull(e_patid)

# Step 3: For each missing patient, find the mode (most frequent value) across imputations
mode_imputed <- combined_imputations %>%
  filter(e_patid %in% missing_patients) %>%
  group_by(e_patid) %>%
  summarise(
    stage_mode = names(which.max(table(stage))),
    .groups = 'drop'
  )

# Step 4: Merge mode-imputed values back to original dataset
dt_final <- dt %>%
  left_join(mode_imputed, by = "e_patid") %>%
  mutate(stage_imp = if_else(stage == "", stage_mode, stage)) %>% 
  select(-stage_mode)
  
saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed.rds"))
  
# add some some variables 
  
pop_nam <- "symp_stageImputed"

dtt <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds")))

dt_final <- dt_final %>% mutate(
  age70plus = if_else(age_index>=70, 1, 0),
  month = ceiling(time2diag/30.5),
  diagnosed_stage = as.integer(stage_imp),
  female = 1 - male, 
  nonwhite = if_else(ethn2=="Non-White", 1, 0)
)
  
saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed.rds"))
