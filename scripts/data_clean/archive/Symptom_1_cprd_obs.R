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

# load symptoms and cancer codelists from Kirsten 
# symp_list <- read_dta(file.path(candetect, "codelists", "symptoms", "symptoms_final_codelist.dta"))
# cancer_code <- read_dta(file.path(candetect, "codelists", "comorbidities", "CPRD_cancer_hx_codelist.dta"))
# folder changed
symp_list <- read_dta(file.path(candetect, "codelists", "Methods","Kirsten&Pav", 
                                "symptoms", "symptoms_final_codelist.dta"))

symp_list <- symp_list %>% 
  select(medcodeid, symptom)

# data table version
# symp_list_dt <- setDT(symp_list)

# match the symptom code from Kirsten to obs data one by one ---------------
# only keep obs with matched symptoms
# save the matched data first and combine them later
# a few mdedcodeids have more than one symptoms
# e.g.
# 1566381000006119 diarrhoea      
# 1566381000006119 rectal_bleed  
# because 1566381000006119 was coded as bloody diarrhoea
# but we treat rectal_bleed and diarrhoea as two symptoms
# so we use a m:m merge strategies

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
      fwrite(rt, file.path(wd, "cprd_obs_symptom_interm", sprintf("cprd_obs_symptom_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)

# Combine all obs files into 1 file ----------------------------------------

folder_path <- file.path(wd, "cprd_obs_symptom_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_symptom_%d_%d.csv", folder_path, 
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

# fwrite(cprd_symp, file = file.path(folder_path, "cprd_obs_symptom.csv"))
saveRDS(cprd_symp, file = file.path(wd, "create_study_pop_interm", "cprd_obs_symptom.rds"))

