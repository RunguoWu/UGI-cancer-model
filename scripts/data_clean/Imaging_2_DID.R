# Clean HES DID data
# Create DID scan file with date, scan_type and site

rm(list = ls())
library(haven)
library(readxl)
library(tidyverse)
library(dtplyr)
library(data.table)
library(bit64)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# read in CT/US and endocopy codes in DID
ct_us <- read_excel(file.path(wd, "code_list", "CT_US_CPRD_DID.xlsx"), sheet = 2)
endo <- read_excel(file.path(wd, "code_list", "endoscopy_CPRD_DID.xlsx"), sheet = 2)

ct_us_endo <- ct_us %>% 
  filter(`KA notes`!= "exclude") %>% 
  select(-`KA notes`) %>% 
  bind_rows(endo) %>% 
  mutate(snomed = as.character(`SCT-ID`)) %>% 
  rename(nicip = `NICIP code`) %>% 
  select(Modality, nicip, snomed)

# Read DID referral -------------------------------------------------------
did_ref <- fread(file.path(link_raw, "e_aurum_hesdid_referral_23_002840_dm.txt")) %>%  
  lazy_dt() %>% 
  mutate(
    # use test request date. if not available, use test request received date
    request_date = if_else(did_date1 != "", did_date1, did_date2),
    request_date = as.Date(request_date, format = "%d/%m/%Y"),
    e_patid = as.character.integer64(e_patid),
    submissiondataid = as.character.integer64(submissiondataid)
  ) %>% 
  select(e_patid, submissiondataid, ic_reftype_desc, did_patsource_code, request_date) %>% 
  as.data.frame()

saveRDS(did_ref, file.path(wd, "create_study_pop_interm", "DID_ref.rds"))


# Read DID test -----------------------------------------------------------
did_test <- fread(file.path(link_raw, "e_aurum_hesdid_test_23_002840_dm.txt")) %>%  
  lazy_dt() %>% 
  mutate(
    # use test date. if not available, use report issue date
    test_date = if_else(did_date3 != "", did_date3, did_date4),
    test_date = as.Date(test_date, format = "%d/%m/%Y"),
    e_patid = as.character.integer64(e_patid),
    submissiondataid = as.character.integer64(submissiondataid),
    did_snomedct_code = as.character.integer64(did_snomedct_code),
    did_snomedct_code = if_else(is.na(did_snomedct_code), "",  did_snomedct_code)
  ) %>% 
  select(e_patid, submissiondataid, did_nicip_code, did_snomedct_code, test_date) %>% 
  as.data.frame()

saveRDS(did_test, file.path(wd, "create_study_pop_interm", "DID_test.rds"))

# combine DID refer and test ----------------------------------------------
did_ref <- readRDS(file.path(wd, "create_study_pop_interm", "DID_ref.rds"))
did_test <- readRDS(file.path(wd, "create_study_pop_interm", "DID_test.rds"))

# extract relevant DID rows
# snomedct code is more complete than nicip.
# missing: snomedct=4318; nicip=2812706
# those with missing snomedct all have nicip

# start with those with snomedct
did_ugi1 <- did_test %>% 
  filter(did_snomedct_code!="") %>% 
  inner_join(ct_us_endo[, c("Modality", "snomed")], 
             join_by(did_snomedct_code==snomed), relationship = "many-to-many") %>% 
  distinct()
# next, those with nicip
did_ugi2 <- did_test %>%
  filter(did_nicip_code !="") %>%
  inner_join(ct_us_endo[, c("Modality", "nicip")],
             join_by(did_nicip_code==nicip), relationship = "many-to-many") %>% 
  distinct()

did_ugi <- did_ugi1 %>% 
  bind_rows(did_ugi2) %>% 
  distinct(e_patid, submissiondataid, test_date, Modality)

did_ugi <- did_ugi %>% 
  left_join(did_ref, by = c("e_patid", "submissiondataid"))

did_ugi <- did_ugi %>% 
  mutate(
    pat_source = case_when(did_patsource_code == 1 ~ "Admitted Patient Care - Inpatient", 
                           did_patsource_code == 2 ~ "Admitted Patient Care - Day case",
                           did_patsource_code == 3 ~ "Out-patient",
                           did_patsource_code == 4 ~ "GP Direct Access",
                           did_patsource_code == 5 ~ "Accident and Emergency Department",
                           did_patsource_code == 6 ~ "Other Health Care Provider",
                           did_patsource_code == 7 ~ "Other",
                           did_patsource_code == 99 ~ "Unkown"
    )
  ) %>% 
  select(-did_patsource_code)

saveRDS(did_ugi, file.path(wd, "create_study_pop_interm", "DID_ugi.rds"))
  
#----
did_ugi <- readRDS(wd, "create_study_pop_interm", "DID_ugi.rds")

rt <- table(did_ugi$ic_reftype_desc, did_ugi$pat_source, useNA = "ifany")
write.csv(rt, file.path(output, "DID_xtab.csv"))

did_ugi$request2test_days <- as.integer(did_ugi$test_date - did_ugi$request_date)





