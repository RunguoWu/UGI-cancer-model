# Clean CPRD Observation level data #
# Extract all CPRD observations and consultations for UGI patients in CPRD and NCRAS

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

# Load first UGI cancer records from NCRAS
ugi_ncras <- readRDS(file = file.path(wd, "create_study_pop_interm", "first_ugi_cancer_ncras.rds"))

# Load CPRD cancer data for exclusion
ugi_cprd <- readRDS(file = file.path(wd, "create_study_pop_interm", "ugi_cancer_dates_cprd.rds"))

ugi <- union(ugi_cprd$e_patid, ugi_ncras$e_patid)

rm(ugi_ncras, ugi_cprd)
# Extract obs and cons from all CPRD obs files ----------------------------

cl <- makePSOCKcluster(16)
clusterExport(cl, c("cprd_raw"))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(dtplyr))
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(bit64))
registerDoSNOW(cl)

cprd_pat <- foreach(i = 0:43) %dopar% {
  
  obs_folder <- (file.path(cprd_raw, "observation_txt"))
  cons_folder <- (file.path(cprd_raw, "consultation_txt"))
  
  for (j in 1:23) {
    
    # Observations
    file_name <- if (j <= 9) 
      sprintf("e_aurum_patlist%d_extract_observation_00%d.txt", i, j) else 
        sprintf("e_aurum_patlist%d_extract_observation_0%d.txt", i, j)
    
    if(file.exists(file.path(obs_folder, file_name))){
      dt <- fread(file.path(obs_folder, file_name)) %>% # data table is faster and it can be used as data frame
        lazy_dt() %>% 
        mutate(
          date = if_else(obsdate != "", obsdate, enterdate),
          date = as.Date(date, format = "%d/%m/%Y"),
          medcodeid = as.character.integer64(medcodeid),
          e_patid = as.character.integer64(e_patid)
        ) %>%
        select(e_patid, date, medcodeid) %>%
        filter(e_patid %in% ugi) %>% 
        as.data.frame()
      
      fwrite(dt, file.path(wd, "cprd_all_obs_cons_cancer_interm", sprintf("cprd_obs_allCancerID_%d_%d.csv", i, j)))
    }
    
    # Consultations
    file_name2 <- if (j <= 9) 
      sprintf("e_aurum_patlist%d_extract_consultation_00%d.txt", i, j) else 
        sprintf("e_aurum_patlist%d_extract_consultation_0%d.txt", i, j) 
    
    if(file.exists(file.path(cons_folder, file_name2))){
      dt2 <- fread(file.path(cons_folder, file_name2)) %>% # data table is faster and it can be used as data frame
        lazy_dt() %>% 
        mutate(
          date = if_else(consdate != "", consdate, enterdate),
          date = as.Date(date, format = "%d/%m/%Y"),
          e_patid = as.character.integer64(e_patid)
        ) %>%
        select(e_patid, date) %>%
        filter(e_patid %in% ugi) %>% 
        as.data.frame()
      
      fwrite(dt2, file.path(wd, "cprd_all_obs_cons_cancer_interm", sprintf("cprd_cons_allCancerID_%d_%d.csv", i, j)))
    }
    
  }
}
stopCluster(cl)

# Combine all obs files into 1 file----------------------------------------
folder_path <- file.path(wd, "cprd_all_obs_cons_cancer_interm")

## All observations for UGI cancer patients----
# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_allCancerID_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
rt <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

fwrite(rt, file = file.path(folder_path, "cprd_obs_allCancerID.csv"))

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
rt <- rt[date < "2025-01-01" & date > "1900-01-01"]
rt <- rt %>% mutate(e_patid=as.character.integer64(e_patid),
                    medcodeid=as.character.integer64(medcodeid))

saveRDS(rt, file = file.path(wd, "create_study_pop_interm", "cprd_obs_allUGI.rds"))

## All consultations for UGI cancer patients----
# Generate all potential file paths
all_files <- sprintf("%s/cprd_cons_allCancerID_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
rt <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

fwrite(rt, file = file.path(folder_path, "cprd_cons_allCancerID.csv"))

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
rt <- rt[date < "2025-01-01" & date > "1900-01-01"]
rt <- rt %>% mutate(e_patid=as.character.integer64(e_patid))

saveRDS(rt, file = file.path(wd, "create_study_pop_interm", "cprd_cons_allUGI.rds"))




