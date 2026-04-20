# Clean CPRD Observation level data #
# For UGI symptoms

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

t2d <- read_dta(file.path(candetect, "codelists", "CPRD_Aurum_Codelist_Library",
                          "Archives", "archive_Jan25", # Update 2026-03-20
                        "Medical_conditions", "CPRD_Aurum_diabetes_type2.dta"))
t2d <- t2d %>% filter(previous != "history") %>% select(medcodeid, condition) %>% 
  rename(symptom = condition)


t1d <- read_dta(file.path(candetect, "codelists", "CPRD_Aurum_Codelist_Library", 
                          "Archives", "archive_Jan25", # Update 2026-03-20
                        "Medical_conditions", "CPRD_Aurum_diabetes_type1.dta"))
t1d <- t1d %>% filter(previous != "history") %>% select(medcodeid, condition) %>% 
  rename(symptom = condition)

symp_list <- rbind(t2d, t1d) %>% distinct()

# data table version
# symp_list_dt <- setDT(symp_list)

# match the symptom code from Kirsten to obs data one by one ---------------
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
      rt <- match_symp(dt, symp_list)
      fwrite(rt, file.path(wd, "cprd_obs_diabetes_interm", sprintf("cprd_obs_diabetes_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)

# Combine all obs files into 1 file ----------------------------------------

folder_path <- file.path(wd, "cprd_obs_diabetes_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_diabetes_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_diabetes <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
cprd_diabetes <- cprd_diabetes[date < "2025-01-01"]

cprd_diabetes <- cprd_diabetes %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_diabetes <- cprd_diabetes %>% mutate(medcodeid=as.character.integer64(medcodeid))
cprd_diabetes <- cprd_diabetes %>% mutate(consid=as.character.integer64(consid))
cprd_diabetes <- cprd_diabetes %>% mutate(obsid=as.character.integer64(obsid))

# fwrite(cprd_symp, file = file.path(folder_path, "cprd_obs_symptom.csv"))
saveRDS(cprd_diabetes, file = file.path(wd, "create_study_pop_interm", "cprd_obs_diabetes.rds"))

