# Clean CPRD Observation level data #
# For all cancer

rm(list = ls())
library(haven)
library(readxl)
library(parallel)
library(doSNOW)
library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# load cancer codelist from WS1
# outdated code list
# cancer_code <- read_dta(file.path(candetect, "codelists", "comorbidities", "CPRD_cancer_hx_codelist.dta"))

# cancer_code <- read_dta(file.path(candetect, "codelists", 
#                                   "CPRD_Aurum_Codelist_Library/Medical_conditions", 
#                                   "Archive",
#                                   "CPRD_Aurum_allcancers.dta")) %>% 
#   select(medcodeid, site)
# site code: 
# 1 bladder
# 2 bone
# 3 brain
# 4 breast
# 5 cervix
# 6 colorectal
# 7 headneck
# 8 kidney
# 9 leukemia
# 10 liver
# 11 lung
# 12 lymphoma
# 13 melanoma
# 14 mesothelioma
# 15 myeloma
# 16 oral
# 17 other
# 18 ovary
## 19 pancreas
# 20 prostate
# 21 testis
# 22 uppergit
# 23 uterine
# 24 ovary or testis
## 25 biliary
## 26 stomach
## 27 esophagus

# WS1 has update CPRD cancer code. 
cancer_code <- read_excel(file.path(candetect, "codelists", 
                                    "CPRD_Aurum_Codelist_Library/Conditions", 
                                    "CPRD_Aurum_allcancers.xlsx")) %>% 
  select(medcodeid, site)

# match the code from Kirsten to obs data one by one ---------------
# only keep obs with matched code
# save the matched data first and combine them later
# use a m:m merge strategies

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
      rt <- match_symp(dt, cancer_code)
      fwrite(rt, file.path(wd, "cprd_obs_cancer_interm", sprintf("cprd_obs_cancer_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)


# Create CPRD cancer variable ---------------------------------------------
## 1. Combine all obs files into 1 file----
folder_path <- file.path(wd, "cprd_obs_cancer_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_cancer_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_cancer <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

fwrite(cprd_cancer, file = file.path(folder_path, "cprd_obs_cancer.csv"))

## 2. Keep the earliest obs date for each type of cancer----
# cprd_cancer <- fread(file.path(folder_path, "cprd_obs_cancer.csv")) 
# Exclude obseravtion date<= "1900-01-01"
# there are a suspiciously large number of "1900-01-01" and "1899-12-31"
# looks like the they are some system-default value for dates
cprd_cancer <- cprd_cancer[date > "1900-01-01"]

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
cprd_cancer <- cprd_cancer[date < "2025-01-01"]

# Keep only the earliest observation per patient and cancer type
cprd_cancer <- cprd_cancer[order(date)][
  , .SD[1], by = .(e_patid, site)]

# cprd_cancer is a data.table obj, ids are saved as integer64,
# convert them into characters
cprd_cancer$e_patid <- as.character.integer64(cprd_cancer$e_patid)
cprd_cancer$consid <- as.character.integer64(cprd_cancer$consid)
cprd_cancer$obsid <- as.character.integer64(cprd_cancer$obsid)
cprd_cancer$medcodeid <- as.character.integer64(cprd_cancer$medcodeid)

saveRDS(cprd_cancer, file = file.path(wd, "create_study_pop_interm", "cprd_cancer_first_obs.rds"))

# Create CPRD cancer first dates ------------------------------------------
# Load CPRD cancer records
cprd_cancer <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_cancer_first_obs.rds"))

# Load in all patient IDs to check
cprd_pat <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_pat.rds"))
cprd_pat$e_patid <- as.character(cprd_pat$e_patid )

### Check all cancer patients are in CPRD?
# cprd_cancer_id <- unique(cprd_cancer$e_patid)
# cprd_pat_id <- unique(cprd_pat$e_patid)
# all(cprd_cancer_id%in%cprd_pat_id)
# TRUE
####

cprd_cancer2 <- cprd_cancer %>% 
  left_join(cprd_pat[, c("e_patid", "dob_imp")]) %>% 
  # filter out those cancer record dates are earlier than birth dates
  filter(date>=dob_imp) %>% 
  # rename to cancerdate
  rename(cancerdate = date)

## create UGI cancer and all cancer first dates----
ugi_dates <- cprd_cancer2 %>% 
  # filter(site %in% c(19, 25, 26, 27)) %>%
  filter(site %in% c("pancreas", "biliary", "stomach", "esophagus")) %>%
  group_by(e_patid) %>%
  summarise(ugi_cancerdate = min(cancerdate, na.rm = TRUE))

all_dates <- cprd_cancer2 %>% 
  group_by(e_patid) %>%
  summarise(all_cancerdate = min(cancerdate, na.rm = TRUE))

saveRDS(ugi_dates, file = file.path(wd, "create_study_pop_interm", "ugi_cancer_dates_cprd.rds"))
saveRDS(all_dates, file = file.path(wd, "create_study_pop_interm", "all_cancer_dates_cprd.rds"))

