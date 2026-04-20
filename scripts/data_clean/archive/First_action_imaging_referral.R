# Create the first action variable and dates
# Combine CPRD imaging, DID and referral

rm(list = ls())
library(tidyverse)
library(data.table)
library(dtplyr)
source("/data/WIPH-CanDetect/HealthEco/route.R")

# Load the incident UGI cancer cohort
ugi_symp <- readRDS(file.path(wd, "ugi_symp.rds"))  

# Referral ----------------------------------------------------------------
# cprd_ref <- readRDS(file.path(he, "Jojo", "Data_Output", "referral_imputed_optimized.rds")) %>% 
#   filter(e_patid %in% ugi_symp$e_patid)

# CPRD referrals with UGI related symptoms records around
cprd_ref_ugi_symp <- readRDS(file.path(wd, "cprd_ref_ugi_symp.rds"))

cprd_ref_clean <- cprd_ref_ugi_symp %>% 
  mutate(refurgencyid = case_when(
    refurgencyid == 1 ~ "2 Week Wait",
    refurgencyid == 2 ~ "Urgent",
    refurgencyid == 3 ~ "Soon",
    refurgencyid == 4 ~ "Routine",
    refurgencyid == 5 ~ "Dated",
    TRUE ~ "Unknown"
  )) %>% mutate(
    action_type = case_when(refurgencyid %in% c("2 Week Wait") ~ "2 Week Wait",
                            refurgencyid %in% c("Urgent","Soon") ~ "Urgent referral",
                            refurgencyid %in% c("Routine","Dated") ~ "Routine referral",
                            refurgencyid == "Unknown" ~ "Unknown referral"
    )
  ) %>% 
  mutate(action_date = as.Date(event_date)) %>% 
  mutate(source = "CPRD referral") %>% 
  select(e_patid, action_date, action_type, source)


# HES OP ------------------------------------------------------------------
# HES OP with relevant specialty seeing
hes_op <- readRDS(file = file.path(wd, "hes_op2.rds"))
# try garth's selection
# hes_op <- readRDS(file = file.path(wd, "hes_op3.rds"))

# summary of the interval between request and appointment
# to impute request date when it is not available.
itv <- hes_op%>% 
  mutate(
    action_type = case_when(priority == 1 ~ "Routine referral",
                            priority == 2 ~ "Urgent referral",
                            priority == 3 ~ "2 Week Wait",
                            priority == 9 ~ "Unknown referral"
    )
  ) %>% 
  filter(reqdate >= "2012-04-01" & reqdate <= "2018-12-31") %>% 
  group_by(action_type) %>% 
  summarise(itv = median(as.numeric(apptdate - reqdate), na.rm = TRUE))

op_record <- hes_op%>% 
  mutate(
    action_type = case_when(priority == 1 ~ "Routine referral",
                            priority == 2 ~ "Urgent referral",
                            priority == 3 ~ "2 Week Wait",
                            priority == 9 ~ "Unknown referral"
    ),
    action_date = if_else(!is.na(reqdate), reqdate, apptdate),
    source = if_else(!is.na(reqdate), "HES OP request", "HES OP appointment")
  ) %>% 
  filter(action_date >= "2012-04-01" & action_date <= "2018-12-31") %>% 
  select(e_patid, action_date, action_type, source)

# adjust action_date for those with only appointment date available
op_record <- op_record %>% 
  mutate(
    action_date_adj = case_when(action_type == "2 Week Wait" & source == "HES OP appointment" ~ 
                                  as.numeric(itv[itv$action_type=="2 Week Wait", "itv"]),
                                action_type == "Routine referral" & source == "HES OP appointment" ~ 
                                  as.numeric(itv[itv$action_type=="Routine referral", "itv"]),
                                action_type == "Urgent referral" & source == "HES OP appointment" ~ 
                                  as.numeric(itv[itv$action_type=="Urgent referral", "itv"]),
                                action_type == "Unknown referral" & source == "HES OP appointment" ~ 
                                  as.numeric(itv[itv$action_type=="Unknown referral", "itv"]),
                                TRUE ~ 0
                                )
  ) %>% 
  mutate(
    action_date = action_date - action_date_adj
  ) %>% 
  select(-action_date_adj)


# Pick up the relevant scan record ----------------------------------------
## DID ----
did_ugi <- readRDS(file = file.path(wd, "DID_ugi.rds"))%>% 
  # only look at cancer patients
  filter(e_patid %in% ugi_symp$e_patid) 

# # Only keep those with a GP referrer or source from GP direct access
# did_ugi_old <- did_ugi %>% 
#   filter(grepl("GP", ic_reftype_desc) | grepl("GP", pat_source)) 

# after meeting updated
# keep all referred by a GP
# + unknown referrer but patient source is "GP Direct access"
did_ugi <- did_ugi %>% 
  filter(grepl("GP", ic_reftype_desc) | (ic_reftype_desc == "Not known" & grepl("GP", pat_source))) 
  
# impute missing request date 
# 1 record = "1971-10-22", while the test date is after 2012-04-01
# consider all request earlier than "2000" are not reliable, recode as NA
# actually, only 1
did_ugi$request_date[did_ugi$request_date<"2000-01-01" & !is.na(did_ugi$request_date)] <- NA

# time from request to test
did_ugi$request2test_days <- as.integer(did_ugi$test_date - did_ugi$request_date) 
# only 4 is negative, recode them as NA
did_ugi$request2test_days[did_ugi$request2test_days<0 & !is.na(did_ugi$request2test_days)] <- NA

# find the median waiting time by scan types
rt <- did_ugi %>% 
  group_by(Modality) %>% summarise(median(request2test_days, na.rm=TRUE))
# Modality                      `median(request2test_days, na.rm = TRUE)`
# <chr>                                                             <dbl>
# 1 Computerized axial tomography                                        14
# 2 Diagnostic ultrasonography                                           26
# 3 Endoscopy                                                            23

# impute missing request date = test date - median waiting days
mod_name <- c("Computerized axial tomography", 
              "Diagnostic ultrasonography",
              "Endoscopy")

for (i in mod_name){
  did_ugi$request_date[is.na(did_ugi$request_date) & 
                         did_ugi$Modality==i] <- 
    as.Date(did_ugi$test_date[is.na(did_ugi$request_date) & 
                                did_ugi$Modality==i] -
    as.numeric(rt[rt$Modality==i, "median(request2test_days, na.rm = TRUE)"]))
}

# only keep the essential information
did_ugi <- did_ugi %>% 
  select(e_patid, request_date, Modality) %>% 
  rename(scan_date = request_date, 
         scan_type = Modality
         ) %>% 
  mutate(
    source = "DID"
  ) %>% 
  left_join(ugi_symp[, c("e_patid", "cancerdate")]) %>% 
  # remove those beyond research scope
  filter(scan_date >= "2012-04-01" & scan_date <= cancerdate) %>% 
  select(-cancerdate)
  
## CPRD scan----
cprd_scan <- readRDS(file = file.path(wd, "cprd_obs_scan.rds")) %>% 
  # only look at cancer patients
  filter(e_patid %in% ugi_symp$e_patid) %>% 
  select(-consid, -obsid, -medcodeid, -surgery) %>% 
  mutate(scan_type = case_when(grepl("ct", scan_type) ~ "CT",
                                grepl("endo", scan_type) ~ "Endoscopy",
                                grepl("us", scan_type) ~ "US",
                                scan_type == "unknown" ~ "Unknown"
  )) %>% 
  rename(scan_date = date) %>%
  mutate(scan_date = as.Date(scan_date)) %>% 
  left_join(ugi_symp[, c("e_patid", "cancerdate")]) %>% 
  # remove those beyond research scope
  filter(scan_date >= "2012-04-01" & scan_date <= cancerdate) %>% 
  select(-cancerdate)

# a few scan types are unknown, but that does not matter
# as we accept various kinds of scans
# the key issue is how to address unknown site


# ### an arbituary decision:----
# # if a unknown site scan happened between a recorded symptom and cancer diagnosis, keep it
# # otherwise, drop it
# unknown_site <- cprd_scan %>%
#   filter(site == "unknown") %>%
#   select(e_patid, date, scan_type)
# 
# unknown_site2 <- ugi_symp %>%
#   select(e_patid, cancerdate, sympdate2yr) %>%
#   inner_join(unknown_site) %>%
#   # only 12567 unknown scans related to incident cancer patients
#   filter(date < cancerdate & date > sympdate2yr & !is.na(sympdate2yr)) %>%
#   select(e_patid, date, scan_type)
# 
# # finally, add back
# cprd_scan_ugi <- cprd_scan %>%
#   filter(site != "unknown") %>%
#   select(e_patid, date, scan_type) %>%
#   bind_rows(unknown_site2) %>%
#   rename(scan_date = date) %>%
#   mutate(scan_date = as.Date(scan_date))

# After meeting, we decide not to try to recover unknown site CT and US
# We keep all unknown site endoscopy

cprd_scan_ugi <- cprd_scan %>% 
  filter(!(site == "unknown" & scan_type != "Endoscopy")) %>% 
  mutate(source = "CPRD") %>% 
  select(e_patid, scan_date, scan_type, source)

## Combine DID and CPRD----
all_scan_ugi <- rbind(did_ugi, cprd_scan_ugi)

all_scan_ugi <- all_scan_ugi %>% 
  mutate(scan_type = case_when(grepl("CT|tomography", scan_type) ~ "CT",
                                grepl("Endo", scan_type) ~ "Endoscopy",
                                grepl("US|ultrasono", scan_type) ~ "US",
                                scan_type == "Unknown" ~ "Unknown scan"
  )) %>% 
  rename(
    action_date = scan_date,
    action_type = scan_type
  )

saveRDS(all_scan_ugi, file.path(wd, "imaging_DID_CPRD.rds"))

## Explore the two data sources for scans----

# # look at those excluded for unknown site CT and US
# cprd_scan2 <- cprd_scan %>%
#   mutate(unsure = if_else(site == "unknown" & scan_type != "Endoscopy", 1, 0)
#   )
# # unique ids 
# dtt1 <- cprd_scan2 %>% 
#   distinct(e_patid)
# 
# # unique ids with at least one known site scan
# dtt2 <- cprd_scan2 %>% 
#   filter(unsure==0) %>% 
#   distinct(e_patid)
# 
# # unique ids with only unknown site scan
# unknown_ids <- setdiff(dtt1$e_patid, dtt2$e_patid)
# 
# length(intersect(unknown_ids, unique(did_ugi$e_patid)))


# Combine referral and scans ----------------------------------------------
first_action <- rbind(all_scan_ugi, cprd_ref_clean, op_record)

first_action <- first_action %>% 
  distinct()

first_action <- first_action %>% 
  left_join(ugi_symp[, c("e_patid", "cancerdate")]) %>% 
  filter(action_date <= cancerdate) %>% 
  select(-cancerdate)

# Keep the last obs/cons before first action ------------------------------

# Load all CPRD observations for all cancer patients
# obs_ugi <- readRDS(file = file.path(wd, "cprd_obs_allUGI.rds"))
cons_ugi <- readRDS(file = file.path(wd, "cprd_cons_allUGI.rds"))
# combine the two as any obs/cons
obs_cons_ugi <- 
  # obs_ugi %>% 
  # select(e_patid, date) %>% 
  # bind_rows(cons_ugi) %>% 
  # distinct(e_patid, date) %>% 
  # filter(e_patid %in% ugi_symp$e_patid) %>% 
  # rename(obs_date = date)
  
  # Now only consider consultation
  # obs have too many records very close to action
  cons_ugi %>% 
  filter(e_patid %in% ugi_symp$e_patid) %>% 
  rename(obs_date = date) %>% 
  distinct()

# create a subset of first actions with the latest consultation before it
first_action <- first_action %>% 
  mutate(action_id = 1:n())

first_action2 <- first_action %>% 
  left_join(obs_cons_ugi, relationship = "many-to-many") %>% 
  # exclude cons within 1 month of action,   
  # as they are highly likely to be from one consultation 
  filter(obs_date<action_date - 30 & obs_date>="2012-04-01") %>% 
  group_by(action_id) %>% 
  mutate(last_cons_B4action_date = max(obs_date, na.rm = TRUE)
  ) %>%
  slice(1) %>%
  ungroup() %>% 
  select(-obs_date) %>% 
  mutate(last_cons_B4action_date = as.Date(last_cons_B4action_date))
  
# summary(as.numeric(first_action2$action_date-first_action2$last_cons_B4action_date))

# merge back to add last consultation dates before action
first_action3 <- first_action %>% 
  left_join(first_action2[, c("action_id", "last_cons_B4action_date")])

first_action3$action_id <- NULL

# saveRDS(first_action3, file.path(wd, "first_action.rds"))

# new version, Garth's selection of HES OP
saveRDS(first_action3, file.path(wd, "first_action_garthHESOP.rds"))
