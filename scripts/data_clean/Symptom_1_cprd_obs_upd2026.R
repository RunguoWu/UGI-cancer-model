# Create the new code lists for CPRD symptoms
# Based on WS1 update in 03/2026

rm(list = ls())
library(haven)
library(parallel)
library(doSNOW)
library(tidyverse)
library(dtplyr)
library(readxl)
library(dplyr)
library(purrr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# 2026 WS1 updated UGI symptoms
# Set the path to your folder
folder_path <- "/data/WIPH-CanDetect/codelists/CPRD_Aurum_Codelist_Library/Symptoms/"

# Get all xlsx files starting with "CPRD_Aurum_"
file_list <- list.files(
  path    = folder_path,
  pattern = "^CPRD_Aurum_.*\\.xlsx$",
  full.names = TRUE
)

# Helper to read one file
read_one <- function(f) {
  df <- read_excel(f)
  
  # If "symptom" column is missing, extract it from the file name
  if (!"symptom" %in% names(df)) {
    symptom_value <- gsub("^.*CPRD_Aurum_(.+)\\.xlsx$", "\\1", basename(f))
    df <- df |> mutate(symptom = symptom_value)
  }
  
  df |> select(medcodeid, symptom)
}

# Apply to all files and stack into one data frame
symp_list <- map(file_list, read_one) |> list_rbind()

# Rename some symptoms to align them to the old version

symp_list <- symp_list %>% 
  mutate(symptom = case_when(
    symptom == "epigastric_mass" ~ "upper_abdo_mass", 
    symptom == "heartburn_reflux_symptoms" ~ "heartburn", 
    symptom == "pruritus" ~ "pruritis", 
    symptom == "weightloss" ~ "weight_loss", 
    symptom == "appetiteloss" ~ "appetite_loss", 
    symptom == "haematemsis" ~ "haematemesis", 
    symptom == "nausea_vomiting" ~ "nausea/vomiting", 
    symptom == "back_pain" ~ "painback", 
    symptom == "fatigue" ~ "fatigue/malaise",
    .default = symptom
  ))

saveRDS(symp_list, file.path(wd, "code_list", "symptom_code_list_WS1202603.rds"))

# match the new updated symptom code to obs data one by one ---------------
# so we use a m:m merge strategies to keep all relevant records

# patlist 0-43
# obs 001-023
cl <- makePSOCKcluster(16)
clusterExport(cl, c("cprd_raw"))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(dtplyr))
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(bit64))
registerDoSNOW(cl)

cprd_pat <- foreach(i = 0:43) %dopar% {
  
  setwd(file.path(cprd_raw, "observation_txt"))
  
  for (j in 1:23) {
    
    file_name <- if (j <= 9) 
      sprintf("e_aurum_patlist%d_extract_observation_00%d.txt", i, j) else 
        sprintf("e_aurum_patlist%d_extract_observation_0%d.txt", i, j)
    
    if(file.exists(file_name)){
      dt <- fread(file_name) # data table is faster and it can be used as data frame
      rt <- match_symp(dt, symp_list, from_date = "2012-04-01")
      fwrite(rt, file.path(wd, "cprd_obs_symptom_interm", sprintf("cprd_obs_symptom_upd2026_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)

# Combine all obs files into 1 file ----------------------------------------

folder_path <- file.path(wd, "cprd_obs_symptom_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_symptom_upd2026_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_symp <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
cprd_symp <- cprd_symp[date < "2025-01-01"]

cprd_symp <- cprd_symp %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_symp <- cprd_symp %>% mutate(medcodeid=as.character.integer64(medcodeid))
cprd_symp <- cprd_symp %>% mutate(consid=as.character.integer64(consid))
cprd_symp <- cprd_symp %>% mutate(obsid=as.character.integer64(obsid))

cprd_symp <- cprd_symp %>% 
  mutate(symptom = case_when(
    symptom == "epigastric_mass" ~ "upper_abdo_mass", 
    symptom == "heartburn_reflux_symptoms" ~ "heartburn", 
    symptom == "pruritus" ~ "pruritis", 
    symptom == "weightloss" ~ "weight_loss", 
    symptom == "appetiteloss" ~ "appetite_loss", 
    symptom == "haematemsis" ~ "haematemesis", 
    symptom == "nausea_vomiting" ~ "nausea/vomiting", 
    symptom == "back_pain" ~ "painback", 
    symptom == "fatigue" ~ "fatigue/malaise",
    .default = symptom
  ))

# fwrite(cprd_symp, file = file.path(folder_path, "cprd_obs_symptom.csv"))
saveRDS(cprd_symp, file = file.path(wd, "create_study_pop_interm", "cprd_obs_symptom_upd2026.rds"))







