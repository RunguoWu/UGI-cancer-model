# Revisit CPRD referrals
# Link to observations
# Only keep referral linked to observation which has related symptom recorded

rm(list = ls())
library(tidyverse)
library(dtplyr)
library(data.table)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Load cancer cohort
ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))

# Load CPRD referrals
cprd_ref <- readRDS(file.path(he, "Jojo", "Data_Output", "referral_imputed_optimized.rds")) %>% 

# Only look at cancer cohort's referrals
cprd_ref_ugi <- cprd_ref %>% 
  filter(e_patid %in% ugi_symp$e_patid) 

# Load CPRD observations with related symptoms
cprd_obs_symptom <- readRDS(file = file.path(wd, "cprd_obs_symptom.rds"))

# Load CPRD observations with abnormal Blood test data
hb <- readRDS(file.path(wd, "haemoglobin_clean.rds")) 
plt <- readRDS(file.path(wd, "platelet_clean.rds")) 

# Combine symptom and blood test ------------------------------------------
hb <- hb %>% mutate(
  symptom = case_when(haemoglobin %in% c("low", "abnormal") ~ "hb_low",
                      # most abnormal hb cases are low hb
                      TRUE ~ "hb_not_low")
) %>% 
  filter(symptom == "hb_low") %>% 
  select(e_patid, obsid, date, symptom)

plt <- plt %>% mutate(
  symptom = case_when(platelet %in% c("high", "abnormal") ~ "plt_high",
                      # high platelet cases are more common in abnormal plt
                      TRUE ~ "plt_not_high")
) %>% 
  filter(symptom == "plt_high") %>% 
  select(e_patid, obsid, date, symptom)

cprd_symp <- cprd_obs_symptom %>% 
  select(e_patid, obsid, date, symptom) %>% 
  bind_rows(hb, plt)

# Pick up referrals with obs id linked directly to symptoms ---------------
obsid_list <- unique(cprd_symp$obsid)

cprd_ref_ugi2 <- cprd_ref_ugi %>% 
  filter(obsid %in% obsid_list) %>% 
  mutate(event_date = as.Date(event_date))

# Only 805 referrals

# Extend to referrals with symptoms around --------------------------------
# 14 days around
cprd_symp2 <- cprd_symp %>% 
  select(e_patid, date) %>% 
  filter(date>= "2012-04-01" & date <= "2018-12-31") %>% 
  rename(symp_date = date) 
  
cprd_ref_ugi3 <- cprd_ref_ugi %>% 
  inner_join(cprd_symp2, relationship = "many-to-many") %>% 
  mutate(event_date = as.Date(event_date),
         symp_date = as.Date(symp_date)) %>% 
  # After merge, one id can have multiple symptom dates
  # Only keep rows where event_date within 14 days around symp_date
  filter(abs(as.numeric(event_date - symp_date)) <= 14) %>% 
  distinct(e_patid, obsid, event_date, .keep_all = TRUE) %>% 
  select(-symp_date)  
  # Actually no need to bind cprd_ref_ugi2, as they have been included
  # but keep the code here anyway

saveRDS(cprd_ref_ugi3, file.path(wd, "cprd_ref_ugi_symp.rds"))

















