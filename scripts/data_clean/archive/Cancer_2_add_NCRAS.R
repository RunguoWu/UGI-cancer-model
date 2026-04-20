# Combine cancer record from CPRD and NCRAS
# The CPRD cancer record is created based on cancer codelist from WS1
# The NCRAS cancer record is directly from WS1

rm(list = ls())

library(tidyverse)
library(haven)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# # Create incident UGI cancer data -----------------------------------------
# # Use NCRAS data only, as CPRD cancer data is not reliable
# # Do not account UGI cancer history, will do it later
# 
# # Load NCRAS cancer records
# # ncras <- read_dta(file.path(candetect, "intermediate_data/ncras", "ncras_cohort_cases.dta"))
# # I prefer to use the all cancer data, without exclusion criteria implemented
# # I will use CPRD cancer cancer data to exclude patients with UGI cancer history
# ncras <- read_dta(file.path(candetect, "intermediate_data/ncras", "ncras_allcancers.dta"))
# # ncras_stage <- read_dta(file.path(candetect, "intermediate_data/ncras", "ncras_stage.dta"))
# 
# # site code in ncras
# # 1 gallbladder
# # 2 pancreas
# # 3 oesophagus
# # 4 stomach
# # 5 bladder
# # 6 bone
# # 7 brain
# # 8 breast
# # 9 cervix
# # 10 headneck
# # 11 kidney
# # 12 leukemia
# # 13 liver
# # 14 lowergit
# # 15 lung
# # 16 lymphoma
# # 17 melanoma
# # 18 mesothelioma
# # 19 myeloma
# # 20 oral
# # 21 ovary
# # 22 prostate
# # 23 testis
# # 24 uterine
# # 25 thyroid
# # 26 other primary
# # 27 secondary
# # 28 poorly defined
# 
# # first, create data with only the first record for each UGI cancer type
# ugi_cancer <- ncras %>% 
#   filter(ugit==1) %>% 
#   select(epatid, cancerdate, site) %>% 
#   group_by(epatid, site) %>% 
#   slice_min(order_by = cancerdate, with_ties = FALSE) %>% 
#   ungroup() %>% 
# # then create a column record the first UGI cancer date 
#   group_by(epatid) %>%
#   mutate(first_ugi_cancerdate = min(cancerdate)) %>%
#   ungroup() %>% 
# # recode site
#   mutate(site = as.numeric(site),
#          site = case_when(
#            site == 1 ~ "galb",
#            site == 2 ~ "panc",
#            site == 3 ~ "oeso",
#            site == 4 ~ "stom"
#          ))
# 
# # nearly 200 individuals with multiple UGI cancer record
# # only keep the earliest one
# ugi_cancer <- ugi_cancer %>% 
#   group_by(epatid) %>% 
#   slice_min(order_by = cancerdate, with_ties = FALSE) %>% 
#   ungroup() %>% 
#   # in this case first_ugi_cancerdate = cancerdate
#   # only keep one
#   select(-first_ugi_cancerdate) %>% 
#   rename(e_patid = epatid) %>% 
#   mutate(e_patid = as.character(e_patid))
# 
# saveRDS(ugi_cancer, file = file.path(wd, "first_ugi_cancer_ncras.rds"))


# Create Incident UGI cancer cohort ---------------------------------------

# Directly use cancer cohort data from WS1
# with stage
ugi_ncras <- read_dta(file.path(candetect, "clean_data/cohort", "cases_all_demographics.dta")) %>% 
  rename(e_patid = epatid)

# 25 patients have multisite cancer
# Load previously coded NCRAS data to check which site happened first
# only keep the first site
ugi_ncras2check <- readRDS(file = file.path(wd, "first_ugi_cancer_ncras.rds")) %>% 
  rename(cancerdate_first = cancerdate)

df <- ugi_ncras %>% filter(multisite==1) %>% 
  left_join(ugi_ncras2check) %>% 
# Create primary site and secondary site relying on date order
  rowwise() %>%
  mutate(
    site_2nd = {
      # Get active cancer sites
      active_sites <- c(
        if(!is.na(gb) && gb == 1) "galb",
        if(!is.na(pancreas) && pancreas == 1) "panc", 
        if(!is.na(oesophagus) && oesophagus == 1) "oeso",
        if(!is.na(gastric) && gastric == 1) "stom"
      )
      # Return the one that's not the primary site
      secondary <- active_sites[active_sites != site]
      if(length(secondary) > 0) secondary[1] else NA_character_
    }
  ) %>%
  ungroup() %>% 
  rename(site_1st = site) %>% 
  select(e_patid, site_1st, site_2nd)

print(df, n=25)

# Combine to create the site variable
ugi_ncras <- ugi_ncras %>% 
  mutate(site = case_when(!is.na(gb) & is.na(multisite) ~ "galb",
                          !is.na(pancreas) & is.na(multisite) ~ "panc",
                          !is.na(oesophagus) & is.na(multisite) ~ "oeso",
                          !is.na(gastric) & is.na(multisite) ~ "stom"
  )) %>% 
  left_join(df) %>% 
  mutate(site = if_else(is.na(site), site_1st, site)) %>% 
  select(-site_1st)

# Load CPRD cancer data for exclusion
ugi_cprd <- readRDS(file = file.path(wd, "ugi_cancer_dates_cprd.rds"))
# # Load first UGI cancer records from NCRAS
# ugi_ncras <- readRDS(file = file.path(wd, "first_ugi_cancer_ncras.rds"))

# consider using 2014-04-01 for cancer record
# leave 2-year window for HES DID data to check imaging record
# HES DID starts from 2012-04-01

# History
ugi_his_cprd <- ugi_cprd %>% 
  filter(ugi_cancerdate< "2014-04-01")

ugi_his_ncras <- ugi_ncras %>% 
  filter(cancerdate< "2014-04-01")

# Incident UGI
# Use NCRAS only for reliability
# exclude those with ugi cancer history
ugi_icd <- ugi_ncras %>% 
  filter(cancerdate>= "2014-04-01") %>% 
  filter(!e_patid %in% union(ugi_his_cprd$e_patid, ugi_his_ncras$e_patid)) %>% 
  rename(stage = stage_ugit) %>% 
  select(-gb, -pancreas, -oesophagus, -gastric)

saveRDS(ugi_icd, file.path(wd, "ugi_icd_stage.rds"))

# export summary table----
ugi_icd <- readRDS(file.path(wd, "ugi_icd_stage.rds"))

rt <- table(ugi_icd$site, ugi_icd$stage)
miss_prob <- round(prop.table(rt, 1)*100,1)
rt2 <- with(ugi_icd[ugi_icd$stage !="",], table(site, stage))
prob <- round(prop.table(rt2, 1)*100,1)

rt0 <- table(ugi_icd$site)

rt3 <- cbind(rt0, rt[,1], miss_prob[,1], rt[,2], prob[,1], rt[,3], prob[,2], rt[,4], 
             prob[,3], rt[,5], prob[,4])

colnames(rt3) <- c("No.", "missing", "%", "stage 1", "%", "stage 2", "%", 
                   "stage 3", "%", "stage 4", "%")

write.csv(rt3, file.path(output, "tab_site_stage.csv"))




