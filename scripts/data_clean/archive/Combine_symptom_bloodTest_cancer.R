# Combine symptom data from CPRD with cancer data
# Include blood test data

rm(list = ls())
library(tidyverse)
library(data.table)
library(dtplyr)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Load data ---------------------------------------------------------------

# Symptom data
cprd_symp <- readRDS(file = file.path(wd, "cprd_obs_symptom.rds"))

# Blood test data
hb <- readRDS(file.path(wd, "haemoglobin_clean.rds")) 
plt <- readRDS(file.path(wd, "platelet_clean.rds")) 

# Load CPRD cancer data
ugi_icd_stage <- readRDS(file.path(wd, "ugi_icd_stage.rds"))

# Combine symptom and blood test ------------------------------------------
hb <- hb %>% mutate(
  symptom = case_when(haemoglobin %in% c("low", "abnormal") ~ "hb_low",
                      # most abnormal hb cases are low hb
                      TRUE ~ "hb_not_low")
) %>% 
  filter(symptom == "hb_low") %>% 
  select(e_patid, date, symptom)

plt <- plt %>% mutate(
  symptom = case_when(platelet %in% c("high", "abnormal") ~ "plt_high",
                      # high platelet cases are more common in abnormal plt
                      TRUE ~ "plt_not_high")
) %>% 
  filter(symptom == "plt_high") %>% 
  select(e_patid, date, symptom)

cprd_symp <- cprd_symp %>% 
  select(e_patid, date, symptom) %>% 
  bind_rows(hb, plt)

# Symptom records for those with UGI cancer diagnosis ---------------------
ugi_symp <- ugi_icd_stage %>% 
  lazy_dt() %>% 
  select(e_patid, cancerdate, site, stage) %>% 
  left_join(cprd_symp, by = "e_patid") %>% 
  mutate(
    sympdate = if_else(date < cancerdate, date, NA),
    # create sympdate2yr as the first symptom recorded in 2 years before diagnosis
    sympdate2yr = if_else(date < cancerdate & date >= cancerdate - 365.25*2, date, NA),
    # also try 1 year
    sympdate1yr = if_else(date < cancerdate & date >= cancerdate - 365.25, date, NA)
  ) %>%
  select(-date) %>% 
  group_by(e_patid) %>%
  mutate(sympdate = min(sympdate, na.rm = TRUE),
         sympdate1yr = min(sympdate1yr, na.rm = TRUE),
         sympdate2yr = min(sympdate2yr, na.rm = TRUE)
         ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    sympdate = as.Date(sympdate),
    sympdate1yr = as.Date(sympdate1yr),
    sympdate2yr = as.Date(sympdate2yr),
    time2diag = as.integer(cancerdate - sympdate),
    time2diag1yr = as.integer(cancerdate - sympdate1yr),
    time2diag2yr = as.integer(cancerdate - sympdate2yr)
  ) %>%
  as.data.frame()

saveRDS(ugi_symp, file.path(wd, "ugi_symp.rds"))

