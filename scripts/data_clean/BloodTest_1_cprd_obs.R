# Clean CPRD Observation level data #
# to identify obs for haemoglobin and platelet

rm(list = ls())
library(haven)
library(parallel)
library(readxl)
library(doSNOW)
library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# specific medcodeids for heamoglobin
hb <- read_xlsx(file.path(candetect, "codelists", 
                         "CPRD_Aurum_Codelist_Library", "Lab_Tests",
                         "cprd_aurum_hb.xlsx"))
hb_medcode <- hb$medcodeid

# specific medcodeids for platelet
platelet <- read_xlsx(file.path(candetect, "codelists", 
                               "CPRD_Aurum_Codelist_Library", "Lab_Tests",
                               "cprd_aurum_platelets.xlsx"))
plt_medcode <- platelet$medcodeid


# Pick up obs with hb and plt medcode -------------------------------------

# patlist 0-43
# obs 001-023
cl <- makePSOCKcluster(8)
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
      dt <- fread(file_name) %>% # data table is faster and it can be used as data frame
        lazy_dt() %>% 
        mutate(
          date = if_else(obsdate != "", obsdate, enterdate),
          date = as.Date(date, format = "%d/%m/%Y"),
          medcodeid = as.character.integer64(medcodeid)
        ) %>%
        select(e_patid, obsid, date, medcodeid, value, numunitid) %>% 
        as.data.frame()
     
      dt_hb <- dt %>% 
        lazy_dt() %>% 
        filter(medcodeid %in% hb_medcode) %>% 
        as.data.frame()
      fwrite(dt_hb, file.path(wd, "cprd_obs_bloodTest_interm", sprintf("cprd_obs_hb_%d_%d.csv", i, j)))
      
      dt_plt <- dt %>% 
        lazy_dt() %>% 
        filter(medcodeid %in% plt_medcode) %>% 
        as.data.frame()
      fwrite(dt_plt, file.path(wd, "cprd_obs_bloodTest_interm", sprintf("cprd_obs_plt_%d_%d.csv", i, j)))
    
    } 
  }
}
stopCluster(cl)


# Combine into single file for hb and plt ---------------------------------

## haemoglobin----
folder_path <- file.path(wd, "cprd_obs_bloodTest_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_hb_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_hb <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
# also only keep those after 2012-04-01
cprd_hb <- cprd_hb[date >= "2012-04-01" & date < "2025-01-01"]
cprd_hb <- cprd_hb %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_hb <- cprd_hb %>% mutate(medcodeid=as.character.integer64(medcodeid))
cprd_hb <- cprd_hb %>% mutate(obsid=as.character.integer64(obsid))

saveRDS(cprd_hb, file = file.path(wd, "create_study_pop_interm", "cprd_obs_hb.rds"))

## platelet----
folder_path <- file.path(wd, "cprd_obs_bloodTest_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_plt_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_plt <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
# also only keep those after 2012-04-01
cprd_plt <- cprd_plt[date >= "2012-04-01" & date < "2025-01-01"]
cprd_plt <- cprd_plt %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_plt <- cprd_plt %>% mutate(medcodeid=as.character.integer64(medcodeid))
cprd_plt <- cprd_plt %>% mutate(obsid=as.character.integer64(obsid))

saveRDS(cprd_plt, file = file.path(wd, "create_study_pop_interm", "cprd_obs_plt.rds"))

