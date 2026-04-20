# Create the study population, 
# using 1 yr before diagnosis as the window to define index

rm(list = ls())

library(tidyverse)
library(haven)
library(data.table)
library(dtplyr)
library(mice)
library(miceadds)
library(parallel)
library(matrixStats)
library(skimr)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# Create the population ---------------------------------------------------
# follow the code in Master_create_study_population.R

## Start from the mid of Step 2 in original code ---------------------------
ugi_icd_stage <- readRDS(file.path(wd, "create_study_pop_interm", "ugi_icd_stage.rds"))
# cprd_symp_hb_plt_diab <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab.rds"))
cprd_symp_hb_plt_diab <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab_upd2026.rds"))
cprd_ppi_ugi <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_ppi_ugi.rds"))

##> Add symptoms ----
ugi_symp <- add_symptom(ugi_icd_stage, cprd_symp_hb_plt_diab)

# saveRDS(ugi_symp, file.path(wd, "ugi_symp.rds"))

##> Additionally add all symptoms within 1 month window since the index date ----
# index=2: look at 2 years back from diagnosis
# index=1: 1 year back
# other: to 2012-04-01
# time_window: time window around index day to include symptoms
ugi_symp2 <- add_symptoms_around_index(ugi_symp, cprd_symp_hb_plt_diab, cprd_ppi_ugi, 
                                       index=1, time_window = 28) # 28 days around index day

##> Save incident UGI cancer cohort with symptom records----
# saveRDS(ugi_symp2, file.path(wd, "ugi_symp_1yrIndex.rds"))
saveRDS(ugi_symp2, file.path(wd, "ugi_symp_1yrIndex_upd2026.rds"))

## Then move on to Step 4 --------------------------------------------------
# Step 3 is irrelevant from the index

##> Load data ----
# Load first action data, including scans and referrals
# first_action <- readRDS(file.path(wd, "create_study_pop_interm", "first_action_cons.rds"))
first_action <- readRDS(file.path(wd, "create_study_pop_interm", "first_action_cons_upd2026.rds"))

# Load the incident UGI cancer cohort
# ugi_symp <- readRDS(file.path(wd, "ugi_symp_1yrIndex.rds"))
ugi_symp <- readRDS(file.path(wd, "ugi_symp_1yrIndex_upd2026.rds"))

# Load CPRD patient ids
cprd_pat <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_pat.rds"))

# Load IMD
imd <- readRDS(file.path(wd, "create_study_pop_interm", "imd.rds"))

# Load WS1 data for ethnicity, deprivation and smoking history
demo <- read_dta(file.path(candetect, "clean_data", "cohort", "cases_all_demographics.dta"))
demo <- demo %>% rename(e_patid = epatid) %>% select(e_patid, ethnicity, smokingstatus)

# Load death data
death <- readRDS(file.path(wd, "create_study_pop_interm", "death.rds"))

##> Combine first action, symptoms and the cohort ----

# First, Patients with symptoms recorded
ugi_symp2 <- combine_symptom_action(ugi_symp, first_action, period = "1yr")

##> Add demographic and death data ----
ugi_symp2 <- add_demo_death(ugi_symp2, cprd_pat, imd, demo, death)

##> Category symptoms ----
# symp_list <- read_dta(file.path(candetect, "codelists", "Methods","Kirsten&Pav", 
#                                 "symptoms", "symptoms_final_codelist.dta"))
symp_list <- readRDS(file.path(wd, "code_list", "symptom_code_list_WS1202603.rds"))
symp_list <- c(unique(symp_list$symptom), "hb_low", "plt_high", 
               "new diabetes T1/2", "tr_dyspepsia")

ng12_symp <- c("diarrhoea", "painback", "constipation", "epigastric_pain", 
               "nausea/vomiting", "dysphagia", "haematemesis", "jaundice",
               "dyspepsia", "weight_loss", "heartburn", "upper_abdo_mass", 
               "hb_low", "plt_high", "new diabetes T1/2", "tr_dyspepsia")

red_flag <- c("dysphagia", "haematemesis", "jaundice", "upper_abdo_mass", "tr_dyspepsia")
ext_symp <- setdiff(symp_list, ng12_symp)

ugi_symp2 <- category_symp(ugi_symp2, symp_list, ng12_symp, red_flag)

# Add NICE NG12 further action criteria
ugi_symp2 <- category_symp_ng12(ugi_symp2)

##> Finally save the study population data ----
# saveRDS(ugi_symp2, file.path(wd, "study_pop_symp_1yrIndex.rds"))
saveRDS(ugi_symp2, file.path(wd, "study_pop_symp_1yrIndex_upd2026.rds"))

# Multiple Imputation -----------------------------------------------------

# Load data to impute
# dt <- readRDS(file.path(wd, "study_pop_symp_1yrIndex.rds"))
dt <- readRDS(file.path(wd, "study_pop_symp_1yrIndex_upd2026.rds"))
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

# saveRDS(imp, file.path(wd, "imp.rds"))
# saveRDS(imp_list, file.path(wd, "imp_list.rds"))

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


# Add variables -----------------------------------------------------------

dt_final <- dt_final %>% mutate(
  age70plus = if_else(age_index>=70, 1, 0),
  month = ceiling(time2diag/30.5),
  diagnosed_stage = as.integer(stage_imp),
  female = 1 - male, 
  nonwhite = if_else(ethn2=="Non-White", 1, 0)
)

dt_final <- dt_final %>% 
  mutate(action_type11 = if_else(action_type %in% c("Routine referral", "Unknown referral"), 
                                 "Routine/unknown referral", action_type),
         action_type11 = factor(action_type11, 
                                levels = c("Endoscopy", "US", "CT", "Unknown scan",
                                           "2 Week Wait", "Routine/unknown referral",
                                           "Urgent referral")
         ),
         site = factor(site, levels = c("oeso", "stom", "panc", "galb"))
  ) %>% 
  mutate(
    ng12_red_flag_match = ifelse(
      (ng12_red_flag %in% c("oeso/stom 2ww", "oeso/stom endoscopy") & site == "oeso") |
        (ng12_red_flag %in% c("oeso/stom 2ww", "oeso/stom endoscopy", "stom/galb 2ww") & site == "stom") |
        (ng12_red_flag %in% c("stom/galb 2ww") & site == "galb") |
        (ng12_red_flag %in% c("panc 2ww", "panc CT/USS") & site == "panc"), "Site-matched red flag",
      ifelse(ng12_red_flag !="No red flag", "Site-unmatched red flag", "No red flag"))
  ) %>% 
  mutate(# meeting decided this categorical variable to be used
    ng12_red_flag2 = ifelse(
      grepl("2ww", ng12_red_flag), "2 Week Wait", 
      ifelse(grepl("endoscopy|CT|USS", ng12_red_flag), "Imaging", "No red flag")
    ),
    ng12_red_flag2 = factor(ng12_red_flag2, levels = c("No red flag", "Imaging", "2 Week Wait"))
  ) %>% 
  mutate(
    action_type21 = factor(action_type21, levels = c("Routine referral", 
                                                     "Urgent referral",
                                                     "Imaging",
                                                     "2 Week Wait"
    ))
  ) %>% 
  mutate(# meeting decided this categorical variable to be used
    ng12_red_flag3 = ifelse(
      grepl("panc", ng12_red_flag), "suspected panc", 
      ifelse(grepl("stom", ng12_red_flag), "suspected oeso/stom/galb", "No red flag")
    ),
    ng12_red_flag3 = factor(ng12_red_flag3, levels = c("No red flag", "suspected panc", "suspected oeso/stom/galb"))
  ) %>% 
  mutate(
    ethnicity = if_else(ethnicity=="", NA, ethnicity), 
    ethnicity = factor(ethnicity, levels = c("White", "Asian", "Black", "Mixed", 
                                             "Other")),
    stage = if_else(stage=="", NA, stage)
  )

dt_final <- dt_final %>% 
  mutate(age10 = cut(
    age_index,
    breaks = c(-Inf, 39.999, 49.999, 59.999, 69.999, 79.999, Inf),
    labels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"),
    right = TRUE
  ),
  factor(age10, levels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"))
  ) %>% 
  mutate(age10_new = case_when(age_index < 60 ~ "<60",
                               age_index < 70 & age_index>=60 ~ "60-69", 
                               age_index < 80 & age_index>=70 ~ "70-79", 
                               age_index >=80 ~ ">=80"
  ),
  age10_new = factor(age10_new, levels = c("<60", "60-69", "70-79", ">=80"))
  )

# create a simpler action variable
dt_final <- dt_final %>% 
  mutate(
    action_type22 = if_else(action_type21=="Urgent referral", "Routine referral", 
                            action_type21),
    action_type22 = factor(action_type22, levels = c("Routine referral", "Imaging", "2 Week Wait"))
  ) 

# change site factor order
dt_final <- dt_final %>% 
  mutate(
    site = factor(site, levels = c("panc", "oeso", "stom", "galb"))
  )

# saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed_1yrIndex.rds"))
saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed_1yrIndex_upd2026.rds"))





