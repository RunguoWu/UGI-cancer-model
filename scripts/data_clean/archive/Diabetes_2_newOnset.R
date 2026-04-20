# Create first diabetes date variable 
# Combine T1 and T2 diabetes

rm(list = ls())

library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")

ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))
cprd_obs_diabetes <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_diabetes.rds"))

ugi_diab <- cprd_obs_diabetes %>% 
  filter(e_patid %in% ugi_symp$e_patid)

# Do not distinguish T1 and T2
# Keep whichever happened first
first_diab <- ugi_diab %>% 
  select(e_patid, obsid, date) %>% 
  filter(date > "1900-01-01") %>% 
  group_by(e_patid) %>%
  arrange(date) %>% 
  mutate(
    diab_date = min(date, na.rm = TRUE)
  ) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(-date) %>% 
  rename(date = diab_date) %>%
  mutate(symptom = "new diabetes T1/2") %>% 
  filter(date >= "2012-04-01") # align with blood test and symptom records

saveRDS(first_diab, file.path(wd, "create_study_pop_interm", "first_diab.rds"))
  
