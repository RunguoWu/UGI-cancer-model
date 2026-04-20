# Clean HES OP data to supplement referral 

rm(list = ls())
library(haven)
library(parallel)
library(doSNOW)
library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# Load the incident UGI cancer cohort
# Only include the cancer cohort
# ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))
ugi_symp <- readRDS(file.path(wd, "ugi_symp_upd2026.rds"))
idlist <- ugi_symp$e_patid

# Only look at 2012-2018
# Because cancer registration data stop at 2018

for (i in 2012:2018) {
  
  setwd(file.path(cprd_raw, "linked_data_txt"))
  
  file_name1 <- sprintf("e_aurum_hesop_patient_pathway_%d_23_002840_dm.txt", i) 
  file_name2 <- sprintf("e_aurum_hesop_appointment_%d_23_002840_dm.txt", i) 
  file_name3 <- sprintf("e_aurum_hesop_clinical_%d_23_002840_dm.txt", i) 
  
  dt1 <- fread(file_name1) %>% 
    lazy_dt() %>% 
    mutate(e_patid=as.character.integer64(e_patid),
           attendkey=as.character.integer64(attendkey)
    ) %>% 
    filter(e_patid %in% idlist) %>% 
    mutate(
      perenddate = as.Date(perend, format = "%d/%m/%Y"),
      perstartdate = as.Date(perstart, format = "%d/%m/%Y"),
      subdate = as.Date(subdate, format = "%d/%m/%Y"),
    ) %>% 
    select(e_patid, attendkey, perenddate, perstartdate, subdate) %>% 
    as.data.frame()
  
  dt2 <- fread(file_name2) %>% 
    lazy_dt() %>% 
    mutate(e_patid=as.character.integer64(e_patid),
           attendkey=as.character.integer64(attendkey)
    ) %>% 
    filter(e_patid %in% idlist) %>% 
    mutate(
      apptdate = as.Date(apptdate, format = "%d/%m/%Y"),
      reqdate = as.Date(reqdate, format = "%d/%m/%Y")
    ) %>% 
    select(e_patid, attendkey, apptdate, attended, firstatt, reqdate, priority, refsourc, waiting) %>% 
    as.data.frame()
  
  dt3 <- fread(file_name3) %>% 
    # lazy_dt() %>% 
    mutate(e_patid=as.character.integer64(e_patid),
           attendkey=as.character.integer64(attendkey)
    ) %>% 
    filter(e_patid %in% idlist) %>% 
    select(e_patid, attendkey, diag_01, opertn_01, operstat, tretspef, mainspef) %>% 
    as.data.frame()
  
  dtt <- dt1 %>% 
    full_join(dt2, by = c("e_patid", "attendkey")) %>% 
    full_join(dt3, by = c("e_patid", "attendkey"))
  
  # fwrite(dtt, file.path(wd, "hes_op_interm", sprintf("hesop_%d.csv", i)))
  fwrite(dtt, file.path(wd, "hes_op_interm", sprintf("hesop_%d_upd2026.csv", i)))
}

# Combine 2012-2018
folder_path <- file.path(wd, "hes_op_interm")

# Generate all potential file paths
all_files <- sprintf("%s/hesop_%d_upd2026.csv", folder_path, 2012:2018)
# Read and combine all existing files efficiently
hes_op <- rbindlist(lapply(all_files, fread), use.names = TRUE, fill = TRUE)

hes_op$e_patid <- as.character.integer64(hes_op$e_patid)
hes_op$attendkey <- as.character.integer64(hes_op$attendkey)
hes_op$perenddate <- as.Date(hes_op$perenddate)
hes_op$perstartdate <- as.Date(hes_op$perstartdate)
hes_op$subdate <- as.Date(hes_op$subdate)
hes_op$apptdate <- as.Date(hes_op$apptdate)
hes_op$reqdate <- as.Date(hes_op$reqdate)

# saveRDS(hes_op, file = file.path(wd, "create_study_pop_interm", "hes_op.rds"))
saveRDS(hes_op, file = file.path(wd, "create_study_pop_interm", "hes_op_upd2026.rds"))
