# Identify PPI prescription

rm(list = ls())
library(haven)
library(parallel)
library(doSNOW)
library(readxl)
library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

ppi <- read_xlsx(file.path(candetect, "codelists", "CPRD_Aurum_Codelist_Library", 
                          "Prescriptions", "CPRD_Aurum_ppi.xlsx"))

ppi_list <- ppi %>% select(prodcodeid)

# match the symptom code from Kirsten to obs data one by one ---------------
cl <- makePSOCKcluster(16)
clusterExport(cl, c("cprd_raw"))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(dtplyr))
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(bit64))
registerDoSNOW(cl)

cprd_pat <- foreach(i = 0:43) %dopar% {
  
  setwd(file.path(cprd_raw, "drug_issue_txt"))
  
  for (j in 1:23) {
    
    file_name <- if (j <= 9) 
      sprintf("e_aurum_patlist%d_extract_drugissue_00%d.txt", i, j) else 
        sprintf("e_aurum_patlist%d_extract_drugissue_0%d.txt", i, j)
    
    if(file.exists(file_name)){
      dt <- fread(file_name) # data table is faster and it can be used as data frame
      rt <- match_drug(dt, ppi_list)
      fwrite(rt, file.path(wd, "cprd_drug_interm", sprintf("cprd_drug_ppi_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)

# Combine all obs files into 1 file ----------------------------------------

folder_path <- file.path(wd, "cprd_drug_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_drug_ppi_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_ppi <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
cprd_ppi <- cprd_ppi[date < "2025-01-01"]

cprd_ppi <- cprd_ppi %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_ppi <- cprd_ppi %>% mutate(prodcodeid=as.character.integer64(prodcodeid))

saveRDS(cprd_ppi, file = file.path(wd, "create_study_pop_interm", "cprd_ppi.rds"))


# Link to cancer cohort ---------------------------------------------------
# Load the incident UGI cancer cohort
ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))

cprd_ppi_ugi <- cprd_ppi %>% 
  filter(cprd_ppi$e_patid %in% ugi_symp$e_patid) %>% 
  filter(date >= "2010-04-01") %>% # do not consider very old PPI prescription
  mutate(ppi_date = as.Date(date)) %>% 
  select(e_patid, ppi_date)

saveRDS(cprd_ppi_ugi, file = file.path(wd, "create_study_pop_interm", "cprd_ppi_ugi.rds"))




