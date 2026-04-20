# Step by Step create the study population for Candetect HE

rm(list = ls())

library(tidyverse)
library(haven)
library(data.table)
library(dtplyr)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# Step 1: Incident UGI cancer cohort --------------------------------------

## 1.0 Create my own NCRAS cohort ----

# I prefer to use the all cancer data, without exclusion criteria implemented
# I will use CPRD cancer cancer data to exclude patients with UGI cancer history
# Use the intermediate data from WS1 to create incident UGI cancer from NCRAS 

##> Load data ----
ncras <- read_dta(file.path(candetect, "intermediate_data/ncras", "ncras_allcancers.dta"))

##> Implement the main function ----
first_ugi_cancer_ncras <- ncras_cancer(ncras)

##> Save NCRAS UGI cohort data ----
# saveRDS(first_ugi_cancer_ncras, file = file.path(wd, "create_study_pop_interm", "first_ugi_cancer_ncras.rds"))
# Note:
# Now, I directly use cleaned NCRAS data from WS1,
# Output above is only used to decide the primary cancer site when multisite cancer happened

## 1.1 Create incident cancer cohort using WS1 data ----

##> Load data ----
# Directly use cancer cohort data from WS1
# with stage
ugi_ncras <- read_dta(file.path(candetect, "clean_data/cohort", "cases_all_demographics.dta")) %>% 
  rename(e_patid = epatid)

# Load pre-prepared CPRD cancer data for exclusion
ugi_cprd <- readRDS(file.path(wd, "create_study_pop_interm", "ugi_cancer_dates_cprd.rds"))

# Load self-made NCRAS data for deciding the primary site of multisite cancer   
first_ugi_cancer_ncras <- readRDS(file.path(wd, "create_study_pop_interm", "first_ugi_cancer_ncras.rds"))

##> Implement the main function ----
ugi_icd <- create_cancer_cohort(ugi_ncras, ugi_cprd, first_ugi_cancer_ncras)

##> Save incident UGI cancer cohort ----
saveRDS(ugi_icd, file.path(wd, "create_study_pop_interm", "ugi_icd_stage.rds"))

# Step 2: Add symptoms to cancer cohort -----------------------------------

##> Load data ----
# Symptom data
# cprd_symp <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_symptom.rds"))
# Update symptom code list from WS1 20260323
cprd_symp <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_symptom_upd2026.rds"))

# Load CPRD observations with related symptoms
hb <- readRDS(file.path(wd, "create_study_pop_interm", "haemoglobin_clean.rds")) 
plt <- readRDS(file.path(wd, "create_study_pop_interm", "platelet_clean.rds")) 

# first_diab <- readRDS(file.path(wd, "create_study_pop_interm", "first_diab.rds"))
# Update symptom code list from WS1 20260323
first_diab <- readRDS(file.path(wd, "create_study_pop_interm", "first_diab_upd2026.rds"))

cprd_ppi_ugi <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_ppi_ugi.rds"))

# Incident UGI cancer data
ugi_icd_stage <- readRDS(file.path(wd, "create_study_pop_interm", "ugi_icd_stage.rds"))

##> Combine symptoms and blood test results----
cprd_symp_hb_plt <- combine_symptom_bt(cprd_symp, hb, plt)
# the intermediate has a huge size, better to save

##> Save a hugh intermiate data----
# saveRDS(cprd_symp_hb_plt, file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt.rds"))

##> add first diabetes ----
# cprd_symp_hb_plt <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt.rds"))

cprd_symp_hb_plt_diab <- bind_rows(cprd_symp_hb_plt, first_diab)

# saveRDS(cprd_symp_hb_plt_diab, file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab.rds"))
# After update WS1 symptom and diabetes code list on 2026-03-23
saveRDS(cprd_symp_hb_plt_diab, file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab_upd2026.rds"))

##> Add symptoms ----
cprd_symp_hb_plt_diab <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab_upd2026.rds"))

ugi_symp <- add_symptom(ugi_icd_stage, cprd_symp_hb_plt_diab)

##> Additionally add all symptoms within 1 month window since the index date ----
# index=2: look at 2 years back from diagnosis
# index=1: 1 year back
# other: to 2012-04-01
# time_window: time window around index day to include symptoms
ugi_symp2 <- add_symptoms_around_index(ugi_symp, cprd_symp_hb_plt_diab, cprd_ppi_ugi, 
                                       index=2, time_window = 28) # 28 days around index day

##> Save incident UGI cancer cohort with symptom records----
# saveRDS(ugi_symp2, file.path(wd, "ugi_symp.rds"))
# After update WS1 symptom and diabetes code list on 2026-03-23
saveRDS(ugi_symp2, file.path(wd, "ugi_symp_upd2026.rds"))

# Step 3: Identify first action -------------------------------------------

##> Load the cohort data prepared in Step 2 ----
# ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))
# Update WS1 symptom and diabetes code list on 2026-03-23
ugi_symp <- readRDS(file.path(wd, "ugi_symp_upd2026.rds"))

## 3.1 CPRD referrals ----

##> Load data ----
# Load extracted CPRD referrals records
cprd_ref <- readRDS(file.path(he, "Jojo", "Data_Output", "referral_imputed_optimized.rds")) 
  # Only look at cancer cohort's referrals
cprd_ref_ugi <- cprd_ref %>% 
  filter(e_patid %in% ugi_symp$e_patid) 

# Load the symptom records to identify relevant referrals
# cprd_symp_hb_plt_diab <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab.rds"))
# Update WS1 symptom and diabetes code list on 2026-03-23
cprd_symp_hb_plt_diab <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_symp_hb_plt_diab_upd2026.rds"))

##> Filter referrals with symptoms ----
# Only keep rows where event_date within 14 days around symp_date
# update 2026-03-23, keep referrals with symptoms 28 days before 
# any referrals with symptoms 7 days after in case input lag from GP
# revised the function in fn_clean.R
cprd_ref_ugi_symp <- cprd_referral_clean(cprd_ref_ugi, cprd_symp_hb_plt_diab, time_window = 28)

##> Save UGI-symptom-related CPRD referrals----
# saveRDS(cprd_ref_ugi_symp, file.path(wd, "create_study_pop_interm", "cprd_ref_ugi_symp.rds"))
# Update WS1 symptom and diabetes code list on 2026-03-23
saveRDS(cprd_ref_ugi_symp, file.path(wd, "create_study_pop_interm", "cprd_ref_ugi_symp_upd2026.rds"))

## 3.2 HES OP for referral ----

##> Load data ----
# extracted HES OP 
# hes_op <- readRDS(file = file.path(wd, "create_study_pop_interm", "hes_op.rds"))
hes_op <- readRDS(file = file.path(wd, "create_study_pop_interm", "hes_op_upd2026.rds")) # actually the same as hes_op.rds, but with additional variables.

# Keep relevant HES OP records
# export mainspef and tretspef value for relevant code identification
# The two variables are for consultant specialty
# rt <- table(hes_op$mainspef, useNA = "ifany")
# write.csv(rt, file.path(output, "HES_OP_mainspef.csv"))
# 
# rt <- table(hes_op$tretspef, useNA = "ifany")
# write.csv(rt, file.path(output, "HES_OP_tretspef.csv"))

# Inclusion code
# refsourc: referral source
# refsourc_code <- c(2, 3, 12)
# mainspef_code <- c("100", "170", "300", "301", "315", "326", "370", "800")
# tretspef_code <- c("100", "104", "105", "106", "170", "173", "300", "301", "306",
#                    "315", "370", "371", "503", "800", "811", "812")

# Garth's: use this
# refsourc_code <- c(2, 3, 12)
# mainspef_code <- c("100", "300", "301", "326")
# tretspef_code <- c("100", "105", "106", "173", "300", "301", "306",
#                    "371", "811", "812")
# write(refsourc_code, file.path(wd, "code_list", "refsourc_code.txt"))
# writeLines(mainspef_code, file.path(wd, "code_list", "mainspef_code.txt"))
# writeLines(tretspef_code, file.path(wd, "code_list", "tretspef_code.txt"))

# Load referral sources codes related to GP
refsourc_code <- scan(file.path(wd, "code_list", "refsourc_code.txt"))
# Load consultant specialty codes related to UGI
mainspef_code <- scan(file.path(wd, "code_list", "mainspef_code.txt"), what = "")
tretspef_code <- scan(file.path(wd, "code_list", "tretspef_code.txt"), what = "")

##> Filter OP records with GP referrals and related specialty ----
op_record <- hes_op_clean(hes_op, refsourc_code, mainspef_code, tretspef_code)

## 3.3 All scans ----

##> Load data ----
# Extracted DID 
did_ugi <- readRDS(file = file.path(wd, "create_study_pop_interm", "DID_ugi.rds"))%>% 
  # only look at the cancer patient cohort
  filter(e_patid %in% ugi_symp$e_patid) 

# Load extracted CPRD scans
cprd_scan <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_scan.rds"))

##> Clean DID data ----
did_ugi2 <- did_clean(did_ugi)

##> Clean CPRD scan data ----
cprd_scan_ugi <- cprd_scan_clean(cprd_scan)

## 3.4 Combine all ----

##> Combine all to create first action data ----
first_action <- create_first_action(did_ugi2, cprd_scan_ugi, cprd_ref_ugi_symp, op_record)

##> Complement using GP consultations for patients with missing symptoms ----
# Add GP consultations only for those without any symptom records

###> Load all CPRD GP consultations for all cancer patients ----
cons_ugi <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_cons_allUGI.rds"))

###> Adding GP consultations into first action ----
first_action_cons <- add_gp_cons(first_action, cons_ugi, gap_days = 28)

##> Save first action data ----
# saveRDS(first_action_cons, file.path(wd, "create_study_pop_interm", "first_action_cons.rds"))
# Update WS1 symptom and diabetes code list on 2026-03-23
saveRDS(first_action_cons, file.path(wd, "create_study_pop_interm", "first_action_cons_upd2026.rds"))

# Step 4: Combine cancer, symptoms, and first actions ---------------------

##> Load data ----
# # Load first action data, including scans and referrals
# first_action <- readRDS(file.path(wd, "create_study_pop_interm", "first_action_cons.rds"))
# # Load the incident UGI cancer cohort
# ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))

# after update WS1 code list on 2026-03-23
first_action <- readRDS(file.path(wd, "create_study_pop_interm", "first_action_cons_upd2026.rds"))
ugi_symp <- readRDS(file.path(wd, "ugi_symp_upd2026.rds"))

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
ugi_symp2 <- combine_symptom_action(ugi_symp, first_action)

# Second, Patients without symptoms recorded
## use_3cons: TRUE = use three consecutive GP consultation as a alert symptom
## FALSE = use the last consultation date before action date as index date 
## cons_interval: define the least intervals between two valid consultations
# ugi_symp3 <- create_index_nosymptom(ugi_symp, first_action, 
#                                     use_3cons=TRUE, cons_interval = 28)
# 
# ugi_symp4 <- create_index_nosymptom(ugi_symp, first_action, 
#                                     use_3cons=FALSE)
# 
# # Third, combine the two and add site and stage
# ugi_symp23 <- rbind(ugi_symp2, ugi_symp3)
# ugi_symp24 <- rbind(ugi_symp2, ugi_symp4)

##> Add demographic and death data ----
ugi_symp2 <- add_demo_death(ugi_symp2, cprd_pat, imd, demo, death)
# ugi_symp3 <- add_demo_death(ugi_symp3, cprd_pat, imd, demo, death)
# ugi_symp4 <- add_demo_death(ugi_symp4, cprd_pat, imd, demo, death)
#   
# ugi_symp23 <- add_demo_death(ugi_symp23, cprd_pat, imd, demo, death)
# ugi_symp24 <- add_demo_death(ugi_symp24, cprd_pat, imd, demo, death)

##> Category symptoms ----
# symp_list <- read_dta(file.path(candetect, "codelists", "Methods","Kirsten&Pav",
#                                 "symptoms", "symptoms_final_codelist.dta"))
# symp_list <- c(unique(symp_list$symptom), "hb_low", "plt_high",
#                "new diabetes T1/2", "tr_dyspepsia")

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
# ugi_symp3 <- category_symp(ugi_symp3, symp_list, ng12_symp, red_flag)
# ugi_symp4 <- category_symp(ugi_symp4, symp_list, ng12_symp, red_flag)
# 
# ugi_symp23 <- category_symp(ugi_symp23, symp_list, ng12_symp, red_flag)
# ugi_symp24 <- category_symp(ugi_symp24, symp_list, ng12_symp, red_flag)

# Add NICE NG12 further action criteria
ugi_symp2 <- category_symp_ng12(ugi_symp2)
# ugi_symp3 <- category_symp_ng12(ugi_symp3)
# ugi_symp4 <- category_symp_ng12(ugi_symp4)
# 
# ugi_symp23 <- category_symp_ng12(ugi_symp23)
# ugi_symp24 <- category_symp_ng12(ugi_symp24)



##> Finally save the study population data ----
# saveRDS(ugi_symp2, file.path(wd, "study_pop_symp.rds"))
saveRDS(ugi_symp2, file.path(wd, "study_pop_symp_upd2026.rds"))

# saveRDS(ugi_symp3, file.path(wd, "study_pop_noSymp_3rdCons.rds"))
# saveRDS(ugi_symp4, file.path(wd, "study_pop_noSymp_ConsB4act.rds"))
# 
# saveRDS(ugi_symp23, file.path(wd, "study_pop_sympPLUSnoSymp_3rdCons.rds"))
# saveRDS(ugi_symp24, file.path(wd, "study_pop_sympPLUSnoSymp_ConsB4act.rds"))


# Step 5 Multiply impute stage --------------------------------------------
library(mice)
library(miceadds)
library(parallel)
library(matrixStats)
library(skimr)

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

# use the mode of imputed values across the 20 imputations as the imputed value ---- 
# Step 1: Combine all imputed datasets with imputation number

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

# add some some variables 
dt_final <- dt_final %>% mutate(
  age70plus = if_else(age_index>=70, 1, 0),
  month = ceiling(time2diag/30.5),
  diagnosed_stage = as.integer(stage_imp),
  female = 1 - male, 
  nonwhite = if_else(ethn2=="Non-White", 1, 0)
)

# saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed.rds"))
saveRDS(dt_final, file.path(wd, "study_pop_symp_stageImputed_upd2026.rds"))

# Recode variables in data for analysis - d4s -----------------------------

# pop_nam <- "symp_stageImputed"
# "symp"
# "noSymp_3rdCons"
# "noSymp_ConsB4act"
# "sympPLUSnoSymp_3rdCons"
# "sympPLUSnoSymp_ConsB4act"

pop_nam <- "symp_stageImputed_upd2026"

d4s <- readRDS(file.path(wd, paste0("study_pop_", pop_nam, ".rds"))) %>% 
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

# Correct age10
d4s <- d4s %>% 
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
d4s <- d4s %>% 
  mutate(
    action_type22 = if_else(action_type21=="Urgent referral", "Routine referral", 
                            action_type21),
    action_type22 = factor(action_type22, levels = c("Routine referral", "Imaging", "2 Week Wait"))
  ) 

# change site factor order
d4s <- d4s %>% 
  mutate(
    site = factor(site, levels = c("panc", "oeso", "stom", "galb"))
  )

saveRDS(d4s, file.path(wd, paste0("study_pop_", pop_nam, ".rds")))
  