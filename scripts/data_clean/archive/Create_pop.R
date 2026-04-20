# Combine symptom, cancer, first action all together
rm(list = ls())
library(tidyverse)
library(data.table)
library(dtplyr)
library(haven)
library(lubridate)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Load first action data, including scans and referrals
# first_action <- readRDS(file.path(wd, "first_action.rds"))

# new version, Garth's selection of HES OP
first_action <- readRDS(file.path(wd, "first_action_garthHESOP.rds"))

# Load the incident UGI cancer cohort
ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))

# Load CPRD patient ids
cprd_pat <- readRDS(file.path(wd, "cprd_pat.rds"))

# Load IMD
imd <- readRDS(file.path(wd, "imd.rds"))

# Load WS1 data for ethnicity, deprivation and smoking history
demo <- read_dta(file.path(candetect, "clean_data", "cohort", "cases_all_demographics.dta"))
demo <- demo %>% rename(e_patid = epatid) %>% select(e_patid, ethnicity, smokingstatus)

# Combine all -------------------------------------------------------------

# first, for those with symptom date, find first action date before diagnosis
# consider 2yr period at this moment
ugi_symp2 <- ugi_symp %>% 
  select(e_patid, cancerdate, sympdate2yr) %>% 
  filter(!is.na(sympdate2yr)) %>% 
  left_join(first_action, by = c("e_patid")) %>% 
  filter(action_date >= sympdate2yr & action_date <= cancerdate) %>% 
  # the filter reduce size from 19564 to 15223
  group_by(e_patid) %>% 
  mutate(first_action_date = min(action_date, na.rm = TRUE)
  ) %>%
  filter(first_action_date == action_date) %>% 
  slice(1) %>%
  ungroup() %>% 
  rename(index_date = sympdate2yr) %>% 
  select(e_patid, index_date, first_action_date, action_type, source, cancerdate)

# second, for those without symptom date, 
# use last consultation date before action date as index date
ugi_symp3 <- ugi_symp %>% 
  select(e_patid, cancerdate, sympdate2yr) %>% 
  filter(is.na(sympdate2yr)) %>% 
  select(-sympdate2yr) %>% 
  left_join(first_action, by = c("e_patid")) %>% 
  filter(action_date >= "2012-04-01" & action_date <= cancerdate) %>% 
  # the filter reduce size from 6884 to 3936
  group_by(e_patid) %>% 
  mutate(first_action_date = min(action_date, na.rm = TRUE)
  ) %>%
  # in this case, need to keep the row with first_action_date and its linked
  # last_cons_B4action_date, so must use the filter below
  filter(first_action_date == action_date) %>% 
  slice(1) %>% 
  ungroup() %>% 
  rename(index_date = last_cons_B4action_date) %>% 
  filter(!is.na(index_date) & index_date >= cancerdate - 365.25*2) %>% 
  # finally left 1708
  # use last consultation b4 action to impute index date
  select(e_patid, index_date, first_action_date, action_type, source, cancerdate)

# third, combine the two and add site and stage
ugi_symp4 <- rbind(ugi_symp2, ugi_symp3)

# # only use those with symptom presentation
# ugi_symp4 <- ugi_symp2

ugi_symp4 <- ugi_symp4 %>% 
  left_join(ugi_symp[, c("e_patid", "site", "stage")])


# Add demographic data ----------------------------------------------------
ugi_symp4 <- ugi_symp4 %>% 
  left_join(subset(cprd_pat, e_patid %in% ugi_symp4$e_patid)) %>% 
  left_join(demo) %>% 
  left_join(subset(imd, e_patid %in% ugi_symp4$e_patid))

ugi_symp5 <- ugi_symp4 %>% mutate(
  time2diag = as.numeric(cancerdate - index_date),
  time2act = as.numeric(first_action_date - index_date),
  action_type2 = case_when(action_type %in% c("CT", "Endoscopy", "US", "Unknown scan") ~ "Imaging",
                           TRUE ~ action_type 
  ),
  action_type3 = case_when(action_type2 == "Imaging" & time2act<=14 ~ "Fast Imaging",
                           action_type2 == "Imaging" & time2act>14 ~ "Slow Imaging",
                           action_type2 == "2 Week Wait" & time2act<=14 ~ "Fast 2 Week Wait",
                           action_type2 == "2 Week Wait" & time2act>14 ~ "Slow 2 Week Wait",
                           (action_type2 == "Urgent referral" | action_type2 == "Unknown referral") & time2act<=14 ~ "Fast urgent referral",
                           action_type2 == "Routine referral" & time2act<=14 ~ "Fast routine referral",
                           (action_type2 == "Urgent referral" | action_type2 == "Unknown referral") & time2act>14 ~ "Slow urgent referral",
                           action_type2 == "Routine referral" & time2act>14 ~ "Slow routine referral",
  ),
  action_type3 = factor(action_type3, levels = c("Fast 2 Week Wait",
                                                 "Fast urgent referral", 
                                                 "Fast routine referral",
                                                 "Fast Imaging",
                                                 "Slow 2 Week Wait",
                                                 "Slow urgent referral",
                                                 "Slow routine referral",
                                                 "Slow Imaging"
  )),
  action_type4 = case_when(action_type2 == "Imaging" & time2act<=14 ~ "Imaging in 2w",
                            action_type2 == "Imaging" & time2act>14 & time2act<=56 ~ "Imaging in 8w",
                           action_type2 == "Imaging" & time2act>56 & time2act<= 182 ~ "Imaging in 26w",
                           action_type2 == "Imaging" & time2act>182 ~ "Imaging after 26w",
                           
                           action_type2 == "2 Week Wait" & time2act<=14 ~ "2 Week Wait in 2w",
                           action_type2 == "2 Week Wait" & time2act>14 & time2act<=56 ~ "2 Week Wait in 8w",
                           action_type2 == "2 Week Wait" & time2act>56 & time2act<=182 ~ "2 Week Wait in 26w",
                           action_type2 == "2 Week Wait" & time2act>182 ~ "2 Week Wait after 26w",
                           
                           action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act<=14 ~ "Non-2ww referral in 2w",
                           action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>14 & time2act<=56  ~ "Non-2ww referral in 8w",
                           action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>56 & time2act<=182  ~ "Non-2ww referral in 26w",
                           action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>182  ~ "Non-2ww referral after 26w",

  ),
  action_type4 = factor(action_type4, levels = c("2 Week Wait in 2w",
                                                  "Imaging in 2w",
                                                  "Non-2ww referral in 2w",
                                                  "2 Week Wait in 8w",
                                                  "Imaging in 8w",
                                                  "Non-2ww referral in 8w",
                                                  "2 Week Wait in 26w",
                                                  "Imaging in 26w",
                                                  "Non-2ww referral in 26w",
                                                  "2 Week Wait after 26w",
                                                  "Imaging after 26w",
                                                  "Non-2ww referral after 26w"
  ))
)  

ugi_symp5 <- ugi_symp5 %>% 
  mutate(age_index = time_length(interval(dob_imp, index_date), "years")) %>% 
  mutate(
    ethn2 = case_when(ethnicity == "White" ~ "White",
                      TRUE ~ "Non-White"),
    ethn2 = factor(ethn2, levels = c("White", "Non-White")),
    age10_cent60 = (age_index - 60)/10) %>% 
  mutate(age10 = cut(
    age_index,
    breaks = c(-Inf, 39, 49, 59, 69, 79, Inf),
    labels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"),
    right = TRUE
  ),
  factor(age10, levels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"))
  ) %>% 
  mutate(imd5_origin = case_when(e2019_imd_10 ==1L | e2019_imd_10 ==2L ~ "1", 
                                 e2019_imd_10 ==3L | e2019_imd_10 ==4L ~ "2",
                                 e2019_imd_10 ==5L | e2019_imd_10 ==6L ~ "3",
                                 e2019_imd_10 ==7L | e2019_imd_10 ==8L ~ "4",
                                 e2019_imd_10 ==9L | e2019_imd_10 ==10L ~ "5"
  )) %>% 
  mutate(
    imd5_imp2 = {
    tbl <- table(imd5_imp, useNA = "no")
    if_else(is.na(imd5_imp), names(tbl)[which.max(tbl)], imd5_imp)
    }
  )

# Add death ---------------------------------------------------------------
death <- readRDS(file.path(wd, "death.rds"))

ugi_symp6 <- ugi_symp5 %>% left_join(death)

saveRDS(ugi_symp6, file.path(wd, "pop_ugi_symp_action_20250903.rds"))


