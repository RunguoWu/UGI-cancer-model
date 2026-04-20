# Clean CPRD Observation level data #
# Imaging scan includes abdo CT & USS and gastroesophageal endoscopy

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

# Load code lists from WS1
endo_code <- read_dta(file.path(candetect, "codelists",  "CPRD_Aurum_Codelist_Library",
                           "Medical_conditions", "CPRD_Aurum_procedures.dta")) %>% 
  filter(condition=="endoscopy") 
# create site column
endo_code <- endo_code %>% 
  mutate(
    site = case_when(grepl("oeso|esoph|reflux|swallow|panendo", term, ignore.case = TRUE) ~ "oesophageal",  
                     grepl("stomach|gastroscopy|gastroduoden", term, ignore.case = TRUE) ~ "gastric",  
                     grepl("upper gi", term, ignore.case = TRUE) ~ "upper GI",
                     grepl("upper gastrointestinal", term, ignore.case = TRUE) ~ "upper GI",
                     grepl("gastrointestinal|digestive", term, ignore.case = TRUE) ~ "GI",
                     grepl("duoden", term, ignore.case = TRUE) ~ "duodenum",  
                     grepl("ileo|ileum", term, ignore.case = TRUE) ~ "ileum",  
                     grepl("jeju", term, ignore.case = TRUE) ~ "jejunum",  
                     grepl("small intest", term, ignore.case = TRUE) ~ "small intestine",
                     .default = "unknown"
                     )
  ) %>% 
  select(medcodeid, term, site) 
# write.csv(endo_code, file.path(output, "endoscopy_code.csv"))

endo_code <- endo_code %>% 
  mutate(scan_type = "endoscopy") %>% 
  select(medcodeid, scan_type, site)

# Load code list from WS1
ct_us_code <- read_excel(file.path(wd, "code_list", "CT_US_CPRD_DID.xlsx"), sheet = 1)
# will use site and scan_type to identify the code of interest
ct_us_code <- ct_us_code %>% 
  select(medcodeid, scan_type, site, surgery) %>% 
  distinct(medcodeid, .keep_all = TRUE)
  
# identify overlap. they are all endoscopic US
dup <- merge(endo_code, ct_us_code, by = "medcodeid")
dup <- dup %>% 
  mutate(scan_type = "endoscopic us", 
         site = if_else(site.y %in% c("esophageal", "gastric"), site.y, site.x)
         ) %>% 
  select(medcodeid, scan_type, site)

# combine all. keep unique medical code
endo_ct_us <- bind_rows(endo_code[!endo_code$medcodeid %in% dup$medcodeid, ], 
                        ct_us_code[!ct_us_code$medcodeid %in% dup$medcodeid, ],
                        dup
                        )


# Match the medical code in endo_ct_us code list --------------------------
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
      rt <- match_symp(dt, endo_ct_us, from_date = "2012-04-01")
      fwrite(rt, file.path(wd, "cprd_obs_scan_interm", sprintf("cprd_obs_scan_%d_%d.csv", i, j)))
    } 
  }
}
stopCluster(cl)

## combine 899 file to 1----
folder_path <- file.path(wd, "cprd_obs_scan_interm")

# Generate all potential file paths
all_files <- sprintf("%s/cprd_obs_scan_%d_%d.csv", folder_path, 
                     rep(0:43, each = 23), rep(1:23, times = 44))
# Filter only files that exist
existing_files <- all_files[file.exists(all_files)]
# Read and combine all existing files efficiently
cprd_scan <- rbindlist(lapply(existing_files, fread), use.names = TRUE, fill = TRUE)

# also there are some dates like 9999-12-31
# only keep those <2025-01-01
# as the data was downloaded before 2025
cprd_scan <- cprd_scan[date < "2025-01-01"]
cprd_scan <- cprd_scan %>% mutate(e_patid=as.character.integer64(e_patid))
cprd_scan <- cprd_scan %>% mutate(medcodeid=as.character.integer64(medcodeid))

cprd_scan <- cprd_scan %>% mutate(site = if_else(site == "esophageal", "oesophageal", site))

# thought
# unknown type of scans are all included, as long as their sites are relevant
# for unknown site, if it follows a related symptom record, include;
# if it does not, include those just before diagnosis

saveRDS(cprd_scan, file = file.path(wd, "create_study_pop_interm", "cprd_obs_scan.rds"))

#----
cprd_scan <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_scan.rds"))

cprd_scan <- cprd_scan %>% 
  mutate(scan_type2 = case_when(grepl("ct", scan_type) ~ "CT",
                                grepl("endo", scan_type) ~ "Endoscopy",
                                grepl("us", scan_type) ~ "US",
                                scan_type == "unknown" ~ "Unknown"
                                ))


rt <- table(cprd_scan$site, cprd_scan$scan_type2,useNA = "ifany")
write.csv(rt, file.path(output, "CPRD_scan_xtab.csv"))




