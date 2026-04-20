# Functions for cleaning Data and create the study population


# CPRD or linkage data extraction -----------------------------------------

# match the symptom code list from Kirsten to CPRD obs level data
match_symp <- function(dt, symp_dic, from_date="1900-01-01"){
  # library(dtplyr)
  rt <- dt %>% 
    lazy_dt() %>% 
    mutate(
      date = if_else(obsdate != "", obsdate, enterdate),
      date = as.Date(date, format = "%d/%m/%Y"),
      medcodeid = as.character.integer64(medcodeid)
    ) %>%
    select(e_patid, consid, obsid, date, medcodeid) %>% 
    filter(date>=from_date) %>% 
    # only keep obs that have matched symptoms
    # and use m:m strategy to keep all symptoms corresponding to one medcodeid
    inner_join(symp_dic, relationship = "many-to-many") %>% 
    as.data.frame()
  
  return(rt)
}

# version dt and symp_dic are data tables
match_symp_dt <- function(dt, symp_dic) {
  # Ensure `symp_dic` is also a data.table
  # setDT(symp_dic)
  
  # Convert column types
  dt[, medcodeid := as.character.integer64(medcodeid)]
  dt[, obsdate := as.IDate(obsdate, format = "%d/%m/%Y")]
  dt[, enterdate := as.IDate(enterdate, format = "%d/%m/%Y")]
  
  # Perform Many-to-Many Join (auto-expanded in `data.table`)
  rt <- merge(
    dt[, .(e_patid, consid, obsid, obsdate, enterdate, medcodeid)], 
    symp_dic, 
    by = "medcodeid", 
    all = FALSE  # Inner join (only matched rows)
  )
  
  return(rt)
}

# compare the two functions.
# actually the first one is a bit faster.
# library(microbenchmark)
# resul <- microbenchmark(
#   match_symp_dt(dt, symp_dic_dt),
#   match_symp(dt, symp_dic),
#   times = 10
# )

# match the symptom code list from Kirsten to CPRD obs level data
match_drug <- function(dt, drug_list, from_date="1900-01-01"){
  # library(dtplyr)
  rt <- dt %>% 
    lazy_dt() %>% 
    mutate(
      date = if_else(issuedate != "", issuedate, enterdate),
      date = as.Date(date, format = "%d/%m/%Y"),
      prodcodeid = as.character.integer64(prodcodeid)
    ) %>%
    select(e_patid, date, prodcodeid) %>% 
    filter(date>=from_date) %>% 
    # only keep obs that have matched symptoms
    # and use m:m strategy to keep all symptoms corresponding to one medcodeid
    inner_join(drug_list, relationship = "many-to-many") %>% 
    as.data.frame()
  
  return(rt)
}


# Creating cancer cohort --------------------------------------------------

ncras_cancer <- function(ncras){

  # site code in ncras
  # 1 gallbladder
  # 2 pancreas
  # 3 oesophagus
  # 4 stomach
  # 5 bladder
  # 6 bone
  # 7 brain
  # 8 breast
  # 9 cervix
  # 10 headneck
  # 11 kidney
  # 12 leukemia
  # 13 liver
  # 14 lowergit
  # 15 lung
  # 16 lymphoma
  # 17 melanoma
  # 18 mesothelioma
  # 19 myeloma
  # 20 oral
  # 21 ovary
  # 22 prostate
  # 23 testis
  # 24 uterine
  # 25 thyroid
  # 26 other primary
  # 27 secondary
  # 28 poorly defined

  # first, create data with only the first record for each UGI cancer type
  ugi_cancer <- ncras %>%
    filter(ugit==1) %>%
    select(epatid, cancerdate, site) %>%
    group_by(epatid, site) %>%
    slice_min(order_by = cancerdate, with_ties = FALSE) %>%
    ungroup() %>%
  # then create a column record the first UGI cancer date
    group_by(epatid) %>%
    mutate(first_ugi_cancerdate = min(cancerdate)) %>%
    ungroup() %>%
  # recode site
    mutate(site = as.numeric(site),
           site = case_when(
             site == 1 ~ "galb",
             site == 2 ~ "panc",
             site == 3 ~ "oeso",
             site == 4 ~ "stom"
           ))

  # nearly 200 individuals with multiple UGI cancer record
  # only keep the earliest one
  first_ugi_cancer_ncras <- ugi_cancer %>%
    group_by(epatid) %>%
    slice_min(order_by = cancerdate, with_ties = FALSE) %>%
    ungroup() %>%
    # in this case first_ugi_cancerdate = cancerdate
    # only keep one
    select(-first_ugi_cancerdate) %>%
    rename(e_patid = epatid) %>%
    mutate(e_patid = as.character(e_patid))

  return(first_ugi_cancer_ncras)
  
}


create_cancer_cohort <- function(ugi_ncras, ugi_cprd, first_ugi_cancer_ncras){
  
  # 25 patients have multisite cancer
  # Load previously coded NCRAS data to check which site happened first
  # only keep the first site
  first_ugi_cancer_ncras <- first_ugi_cancer_ncras %>% 
    rename(cancerdate_first = cancerdate)
  
  df <- ugi_ncras %>% filter(multisite==1) %>% 
    left_join(first_ugi_cancer_ncras) %>% 
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
  
  return(ugi_icd)
} 

# Summarise the missingness across sites and stages  
sum_miss_site_stage <- function(ugi_icd){
  
  rt <- table(ugi_icd$site, ugi_icd$stage)
  miss_prob <- round(prop.table(rt, 1)*100,1)
  rt2 <- with(ugi_icd[ugi_icd$stage !="",], table(site, stage))
  prob <- round(prop.table(rt2, 1)*100,1)
  
  rt0 <- table(ugi_icd$site)
  
  rt3 <- cbind(rt0, rt[,1], miss_prob[,1], rt[,2], prob[,1], rt[,3], prob[,2], rt[,4], 
               prob[,3], rt[,5], prob[,4])
  
  colnames(rt3) <- c("No.", "missing", "%", "stage 1", "%", "stage 2", "%", 
                     "stage 3", "%", "stage 4", "%")
  
  return(rt3)
}



# Add symptoms ------------------------------------------------------------
# symptoms include two abnormal blood test results 

combine_symptom_bt <- function(cprd_symp, hb, plt){
  
  hb <- hb %>% mutate(
    symptom = case_when(haemoglobin %in% c("low", "abnormal") ~ "hb_low",
                        # most abnormal hb cases are low hb
                        TRUE ~ "hb_not_low")
  ) %>% 
    filter(symptom == "hb_low") %>% 
    select(e_patid, obsid, date, symptom)
  
  plt <- plt %>% mutate(
    symptom = case_when(platelet %in% c("high", "abnormal") ~ "plt_high",
                        # high platelet cases are more common in abnormal plt
                        TRUE ~ "plt_not_high")
  ) %>% 
    filter(symptom == "plt_high") %>% 
    select(e_patid, obsid, date, symptom)
  
  cprd_symp <- cprd_symp %>% 
    select(e_patid, obsid, date, symptom) %>% 
    bind_rows(hb, plt)
  
  return(cprd_symp)
}

# cprd_symp_hb_plt is the output of the function above
# cprd_symp_hb_plt_diab is cprd_symp_hb_plt adding diabetes
add_symptom <- function(ugi_icd_stage, cprd_symp_hb_plt_diab){
  
  cprd_symp <- cprd_symp_hb_plt_diab %>% 
    select(e_patid, date, symptom) 
  
  # Add symptom records for those with UGI cancer diagnosis
  ugi_symp <- ugi_icd_stage %>% 
    lazy_dt() %>% 
    select(e_patid, cancerdate, site, stage) %>% 
    left_join(cprd_symp, by = "e_patid") %>% 
    mutate(
      sympdate = if_else(date < cancerdate, date, NA),
      # create sympdate2yr as the first symptom recorded in 2 years before diagnosis
      sympdate2yr = if_else(date < cancerdate & date >= cancerdate - 365.25*2, date, NA),
      # also try 1 year
      sympdate1yr = if_else(date < cancerdate & date >= cancerdate - 365.25, date, NA)
    ) %>%
    select(-date) %>% 
    group_by(e_patid) %>%
    mutate(sympdate = min(sympdate, na.rm = TRUE),
           sympdate1yr = min(sympdate1yr, na.rm = TRUE),
           sympdate2yr = min(sympdate2yr, na.rm = TRUE)
    ) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      sympdate = as.Date(sympdate),
      sympdate1yr = as.Date(sympdate1yr),
      sympdate2yr = as.Date(sympdate2yr),
      # time2diag = as.integer(cancerdate - sympdate),
      # time2diag1yr = as.integer(cancerdate - sympdate1yr),
      # time2diag2yr = as.integer(cancerdate - sympdate2yr)
    ) %>%
    as.data.frame()
  
  return(ugi_symp)
}

# Add all symptoms within 1 month window since the index date
add_symptoms_around_index <- function(ugi_symp, cprd_symp_hb_plt_diab, cprd_ppi_ugi,
                                      index=2, time_window = 14){
  
  if(index==2){
    reference_date <- "sympdate2yr"
  } else if(index==1){
    reference_date <- "sympdate1yr"
  } else {
    reference_date <- "sympdate"
  }
  
  cprd_symp <-  cprd_symp_hb_plt_diab %>% 
    lazy_dt() %>%
    left_join(cprd_ppi_ugi, by = "e_patid", relationship = "many-to-many") %>%
    # recode "dyspepsia" to "tr_dyspepsia" - treatment resistent for qualifying symptoms
    # refer to https://doi.org/10.1016/j.canep.2021.101969 Price 2021
    mutate(ppi_days_diff = as.numeric(as.Date(date) - ppi_date),
           is_tr = (symptom == "dyspepsia" | symptom == "belching_wind") & # add "belching_wind" after update WS1 code list in 2026-03
             !is.na(ppi_days_diff) &
             ppi_days_diff >= 8*7 & ppi_days_diff <=365) %>% 
    group_by(e_patid, date, symptom) %>% 
    summarise(
      has_qualifying_ppi = any(is_tr, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(
      symptom = if_else(has_qualifying_ppi, "tr_dyspepsia", symptom)
    ) %>% 
    select(e_patid, date, symptom) %>% 
    as.data.frame()
  
  # Merge symptom data with cohort
  # Consider drug PPI prescription within 1 year to at least 8 weeks before dyspepsia
  # code it as tr_dyspepsia - treatment resistant
  symptoms_filtered <- ugi_symp %>%
    lazy_dt() %>%
    select(-symptom) %>% # the old symptom names are useless, as they are not specific to time windows
    left_join(cprd_symp, by = "e_patid") %>% 
    # Step 2: Filter symptoms within 28-day window of chosen reference date
    # Using -14 to +14 days window (28 days total centered on reference date)
    filter(!is.na(date)) %>%  # Remove patients without symptom records
    mutate(
      ref_date = get(reference_date),  # Get the chosen reference date
      days_diff = as.numeric(difftime(date, ref_date, units = "days"))
    ) %>%
    filter(days_diff >= -time_window & days_diff <= time_window & date < cancerdate) %>%  # 28-day window: ±14 days
    select(e_patid, symptom, date, days_diff) %>% 
    as.data.frame()
  
  symptom_list <- symptoms_filtered %>%
    group_by(e_patid) %>%
    summarise(symptoms_28day = paste(unique(symptom), collapse = ", ")) %>%
    ungroup()
  
  symptom_total_count <- symptoms_filtered %>%
    group_by(e_patid) %>%
    summarise(total_symptoms_28day = n()) %>%
    ungroup()
  
  ugi_symp2 <- ugi_symp %>%
    left_join(symptom_list, by = "e_patid") %>% 
    left_join(symptom_total_count, by = "e_patid")
  
  return(ugi_symp2)
}

# Define referrals --------------------------------------------------------

cprd_referral_clean <- function(cprd_ref_ugi, cprd_symp_hb_plt_diab, time_window = 28){
  
  # Pick up referrals with obs id linked directly to symptoms 
  obsid_list <- unique(cprd_symp_hb_plt_diab$obsid)
  
  # cprd_ref_ugi2 <- cprd_ref_ugi %>% 
  #   filter(obsid %in% obsid_list) %>% 
  #   mutate(event_date = as.Date(event_date))
  
  # Only 805 referrals
  
  # Extend to referrals with symptoms around
  # 14 days around
  # after discussion extended to 28 days before and 7 days after referral
  cprd_symp2 <- cprd_symp_hb_plt_diab %>% 
    select(e_patid, date) %>% 
    filter(date>= "2012-04-01" & date <= "2018-12-31") %>% 
    rename(symp_date = date) 
  
  cprd_ref_ugi3 <- cprd_ref_ugi %>% 
    inner_join(cprd_symp2, relationship = "many-to-many") %>% 
    mutate(event_date = as.Date(event_date),
           symp_date = as.Date(symp_date)) %>% 
    # After merge, one id can have multiple symptom dates
    # filter(abs(as.numeric(event_date - symp_date)) <= time_window) %>% 
    # update 2026-03-23, keep referrals with symptoms 14 days before 
    # any referrals with symptoms 7 days after in case input lag from GP
    filter(as.numeric(event_date - symp_date) <= time_window &
           as.numeric(event_date - symp_date) >= -(time_window/4) 
           ) %>% 
    distinct(e_patid, obsid, event_date, .keep_all = TRUE) %>% 
    select(-symp_date) %>% 
    mutate(has_symptom = "with symptom")
  # Actually no need to bind cprd_ref_ugi2, as they have been included
  # but keep the code here anyway
  # 
  # update 2026-03 keep all 2WW referral
  cprd_ref_ugi4 <- cprd_ref_ugi %>%
    filter(refurgencyid == 1) %>% 
    mutate(has_symptom = "without symptom")

  cprd_ref_ugi5 <- rbind(cprd_ref_ugi3, cprd_ref_ugi4) %>%
    distinct(e_patid, obsid, event_date, e_pracid, refsourceorgid,reftargetorgid, refurgencyid, refservicetypeid, .keep_all = TRUE)

  return(cprd_ref_ugi5)
  
  # return(cprd_ref_ugi3)
}

hes_op_clean <- function(hes_op, refsourc_code, 
                         mainspef_code, tretspef_code, 
                         all_2ww = FALSE){
  
  # HES OP with relevant specialty seeing
  # hes_op <- readRDS(file = file.path(wd, "hes_op2.rds"))
  # try garth's selection
  if (all_2ww) {
    hes_op <- hes_op %>% 
      filter(refsourc %in% refsourc_code) %>% 
      filter((mainspef %in% mainspef_code) | (tretspef %in% tretspef_code) | priority == 3)
  } else {
    hes_op <- hes_op %>% 
      filter(refsourc %in% refsourc_code) %>% 
      filter((mainspef %in% mainspef_code) | (tretspef %in% tretspef_code))
  }

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
    select(e_patid, action_date, action_type, source, refsourc, mainspef, tretspef)
  
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
  
  return(op_record)
}

# don't filter by specialty
# only filter source - GP, date, and non-first attendance
hes_op_clean2 <- function(hes_op 
                         # mainspef_code, tretspef_code, 
                         #all_2ww = FALSE
                         ){
  
    # hes_op <- hes_op %>% 
    #   filter(refsourc %in% refsourc_code) #%>% # from GP only
    #   # filter(firstatt %in% c("1", "3", "X")) %>%  # first attendance only or unknown
    #   # filter(attended %in% c(5,6,9)) # only attended or unknown
      
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
    select(e_patid, action_date, apptdate, firstatt, attended, action_type, source, refsourc, mainspef, tretspef)
  
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
  
  return(op_record)
}





# All scans ---------------------------------------------------------------

did_clean <- function(did_ugi){
  
  if (!exists("ugi_symp", envir = .GlobalEnv)) {
    assign("ugi_symp", readRDS(file.path(wd, "ugi_symp.rds")), envir = .GlobalEnv)
  } 
  
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
  
  return(did_ugi)
}

cprd_scan_clean <- function(cprd_scan){
  
  if (!exists("ugi_symp", envir = .GlobalEnv)) {
    assign("ugi_symp", readRDS(file.path(wd, "ugi_symp.rds")), envir = .GlobalEnv)
  } 
  
  cprd_scan <- cprd_scan %>% 
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
  
  # ### an arbituary decision:
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
  
  return(cprd_scan_ugi)
}


# Create first action data ------------------------------------------------

create_first_action <- function(did_ugi, cprd_scan_ugi, cprd_ref_ugi_symp, op_record){
  
  if (!exists("ugi_symp", envir = .GlobalEnv)) {
    assign("ugi_symp", readRDS(file.path(wd, "ugi_symp.rds")), envir = .GlobalEnv)
  } 
  
  # Combine DID and CPRD scans
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
  
  # recode CPRD referral variables
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
  
  # combine all 
  first_action <- rbind(all_scan_ugi, cprd_ref_clean, op_record)
  
  first_action <- first_action %>% 
    distinct()
  
  first_action <- first_action %>% 
    left_join(ugi_symp[, c("e_patid", "cancerdate")]) %>% 
    filter(action_date <= cancerdate) %>% 
    select(-cancerdate)
  
  return(first_action)
}

# Add GP consultations only for those without any symptom records
add_gp_cons <- function(first_action, cons_ugi, gap_days = 28){
  
  if (!exists("ugi_symp", envir = .GlobalEnv)) {
   assign("ugi_symp", readRDS(file.path(wd, "ugi_symp.rds")), envir = .GlobalEnv)
  } 
  
  # # Load all CPRD observations for all cancer patients
  # # obs_ugi <- readRDS(file = file.path(wd, "cprd_obs_allUGI.rds"))
  # cons_ugi <- readRDS(file = file.path(wd, "cprd_cons_allUGI.rds"))
  
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
    filter(obs_date < action_date - gap_days & obs_date>="2012-04-01") %>% 
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
  
  return(first_action3)
}


# Combine all -------------------------------------------------------------

# For patients with symptoms recorded
# Use the first symptom as the index date
combine_symptom_action <- function(ugi_symp, first_action, period = "2yr"){

  # Validate period argument
  if (!period %in% c("1yr", "2yr")) {
    stop("period must be either '1yr' or '2yr'")
  }
  
  # Determine which symptom date column to use
  sympdate_col <- paste0("sympdate", period)
  
  # find first action date before diagnosis
  # consider 2yr period at this moment
  ugi_symp2 <- ugi_symp %>% 
    select(e_patid, cancerdate, all_of(sympdate_col)) %>% 
    filter(!is.na(.data[[sympdate_col]])) %>% 
    left_join(first_action, by = c("e_patid")) %>% 
    filter(action_date >= .data[[sympdate_col]] & action_date <= cancerdate) %>% 
    # the filter reduce size from 19564 to 15223
    group_by(e_patid) %>% 
    mutate(first_action_date = min(action_date, na.rm = TRUE)) %>%
    filter(first_action_date == action_date) %>% 
    slice(1) %>%
    ungroup() %>% 
    rename(index_date = all_of(sympdate_col)) %>% 
    select(e_patid, index_date, first_action_date, action_type, source, cancerdate)
  
  ugi_symp2 <- ugi_symp2 %>% 
    left_join(ugi_symp[, c("e_patid", "site", "stage", "symptoms_28day", "total_symptoms_28day")])
  
  return(ugi_symp2)
}

# For patients without symptoms recorded
# Use consultation dates to create the index date
create_index_nosymptom <- function(ugi_symp, first_action, 
                                   use_3cons = TRUE, cons_interval = 0){
  
  # For those without symptoms recorded
  # use_3cons: TRUE = use three consecutive GP consultation as a alert symptom
  # FALSE = use the last consultation date before action date as index date 
  # cons_interval: define the least intervals between two valid consultations
  
  if (!use_3cons) {
    # use last consultation date before action date as index date
    ugi_symp3 <- ugi_symp %>% 
      select(e_patid, cancerdate, sympdate2yr) %>% 
      filter(is.na(sympdate2yr)) %>% 
      select(-sympdate2yr) %>% 
      left_join(first_action, by = c("e_patid")) %>% 
      filter(action_date >= "2012-04-01" & action_date <= cancerdate) %>% 
      # the filter reduce size from 6884 to 3936
      group_by(e_patid) %>% 
      mutate(first_action_date = min(action_date, na.rm = TRUE)
      ) %>%
      # in this case, need to keep the row with first_action_date and its linked
      # last_cons_B4action_date, so must use the filter below
      filter(first_action_date == action_date) %>% 
      slice(1) %>% 
      ungroup() %>% 
      rename(index_date = last_cons_B4action_date) %>% 
      filter(!is.na(index_date) & index_date >= cancerdate - 365.25*2) %>% 
      # finally left 1708
      # use last consultation b4 action to impute index date
      select(e_patid, index_date, first_action_date, action_type, source, cancerdate)
    
  } else {
    # Keep those with at least three GP consultations
    # Select the third one as the index date
    # ref: https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(12)70041-4/abstract
    
    # Load and clean UGI patients GP consultation records
    if (!exists("cons_ugi", envir = .GlobalEnv)) {
      assign("cons_ugi", readRDS(file.path(wd, "create_study_pop_interm", "cprd_cons_allUGI.rds")), envir = .GlobalEnv)
    } 
    
    # Exclude consultations that happened within X days after the last consultation
    # X = cons_interval
    # consultations within this interval are likely strongly associated
    
    obs_cons_ugi <- cons_ugi %>%
      filter(e_patid %in% ugi_symp$e_patid) %>% 
      rename(obs_date = date) %>%
      mutate(obs_date = as.Date(obs_date)) %>% 
      arrange(e_patid, obs_date) %>%
      group_by(e_patid) %>%
      mutate(
        days_since_last = as.numeric(obs_date - lag(obs_date)),
        keep = is.na(days_since_last) | days_since_last >= cons_interval
      ) %>%
      filter(keep) %>%
      select(-days_since_last, -keep) %>%
      ungroup() %>% 
      distinct()
    
    ugi_symp3 <- ugi_symp %>% 
      select(e_patid, cancerdate, sympdate2yr) %>% 
      filter(is.na(sympdate2yr)) %>% 
      select(-sympdate2yr) %>% 
      left_join(first_action, by = c("e_patid")) %>% 
      filter(action_date >= cancerdate - 365.25*2 & action_date <= cancerdate) %>% 
      group_by(e_patid) %>% 
      mutate(first_action_date = min(action_date, na.rm = TRUE)
      ) %>%
      filter(first_action_date == action_date) %>% 
      slice(1) %>% 
      ungroup() %>% 
      # Add all GP consultation dates
      # Keep those with at least three GP consultations
      # Select the third one as the index date
      left_join(obs_cons_ugi, relationship = "many-to-many") %>% 
      # Only keep consultations before first action date and within 2 years before diagnosis
      filter(obs_date <= first_action_date & obs_date >= cancerdate - 365.25*2) %>% 
      arrange(e_patid, obs_date) %>% 
      group_by(e_patid) %>% 
      filter(n() >= 3) %>% 
      mutate(third_obs_date = nth(obs_date, 3)) %>%
      filter(obs_date == third_obs_date) %>% 
      slice(1) %>% 
      ungroup() %>% 
      rename(index_date = third_obs_date) %>% 
      # mutate(index_date = as.Date(index_date)) %>% 
      select(e_patid, index_date, first_action_date, action_type, source, cancerdate)
  }
  
  ugi_symp3 <- ugi_symp3 %>% 
    left_join(ugi_symp[, c("e_patid", "site", "stage", "symptoms_28day", "total_symptoms_28day")])
  
  return(ugi_symp3)
} 

add_demo_death <- function(ugi_symp4, cprd_pat, imd, demo, death){
  
  ugi_symp4 <- ugi_symp4 %>% 
    left_join(subset(cprd_pat, e_patid %in% ugi_symp4$e_patid)) %>% 
    left_join(demo) %>% 
    left_join(subset(imd, e_patid %in% ugi_symp4$e_patid))
  
  ugi_symp5 <- ugi_symp4 %>% mutate(
    time2diag = as.numeric(cancerdate - index_date),
    time2act = as.numeric(first_action_date - index_date),
    time_act2diag = as.numeric(cancerdate - first_action_date),
    action_type2 = case_when(action_type %in% c("CT", "Endoscopy", "US", "Unknown scan") ~ "Imaging",
                             TRUE ~ action_type 
    ),
    action_type21 = if_else(action_type2 =="Unknown referral", "Routine referral", action_type2),
    action_type3 = case_when(action_type2 == "Imaging" & time2act<=14 ~ "Fast Imaging",
                             action_type2 == "Imaging" & time2act>14 ~ "Slow Imaging",
                             action_type2 == "2 Week Wait" & time2act<=14 ~ "Fast 2 Week Wait",
                             action_type2 == "2 Week Wait" & time2act>14 ~ "Slow 2 Week Wait",
                             (action_type2 == "Urgent referral" | action_type2 == "Unknown referral") & time2act<=14 ~ "Fast urgent referral",
                             action_type2 == "Routine referral" & time2act<=14 ~ "Fast routine referral",
                             (action_type2 == "Urgent referral" | action_type2 == "Unknown referral") & time2act>14 ~ "Slow urgent referral",
                             action_type2 == "Routine referral" & time2act>14 ~ "Slow routine referral",
    ),
    action_type3 = factor(action_type3, levels = c("Fast 2 Week Wait",
                                                   "Fast urgent referral", 
                                                   "Fast routine referral",
                                                   "Fast Imaging",
                                                   "Slow 2 Week Wait",
                                                   "Slow urgent referral",
                                                   "Slow routine referral",
                                                   "Slow Imaging"
    )),
    action_type4 = case_when(action_type2 == "Imaging" & time2act<=14 ~ "Imaging in 2w",
                             action_type2 == "Imaging" & time2act>14 & time2act<=56 ~ "Imaging in 8w",
                             action_type2 == "Imaging" & time2act>56 & time2act<= 182 ~ "Imaging in 26w",
                             action_type2 == "Imaging" & time2act>182 ~ "Imaging after 26w",
                             
                             action_type2 == "2 Week Wait" & time2act<=14 ~ "2 Week Wait in 2w",
                             action_type2 == "2 Week Wait" & time2act>14 & time2act<=56 ~ "2 Week Wait in 8w",
                             action_type2 == "2 Week Wait" & time2act>56 & time2act<=182 ~ "2 Week Wait in 26w",
                             action_type2 == "2 Week Wait" & time2act>182 ~ "2 Week Wait after 26w",
                             
                             action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act<=14 ~ "Non-2ww referral in 2w",
                             action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>14 & time2act<=56  ~ "Non-2ww referral in 8w",
                             action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>56 & time2act<=182  ~ "Non-2ww referral in 26w",
                             action_type2 %in% c("Urgent referral", "Unknown referral", "Routine referral") & time2act>182  ~ "Non-2ww referral after 26w",
                             
    ),
    action_type4 = factor(action_type4, levels = c("2 Week Wait in 2w",
                                                   "Imaging in 2w",
                                                   "Non-2ww referral in 2w",
                                                   "2 Week Wait in 8w",
                                                   "Imaging in 8w",
                                                   "Non-2ww referral in 8w",
                                                   "2 Week Wait in 26w",
                                                   "Imaging in 26w",
                                                   "Non-2ww referral in 26w",
                                                   "2 Week Wait after 26w",
                                                   "Imaging after 26w",
                                                   "Non-2ww referral after 26w"
    ))
  )  
  
  ugi_symp5 <- ugi_symp5 %>% 
    mutate(age_index = time_length(interval(dob_imp, index_date), "years")) %>% 
    mutate(
      ethn2 = case_when(ethnicity %in% c("White", "") ~ "White",
                        TRUE ~ "Non-White"),
      ethn2 = factor(ethn2, levels = c("White", "Non-White")),
      age10_cent60 = (age_index - 60)/10) %>% 
    # mutate(age10 = cut( 
    #   age_index,
    #   breaks = c(-Inf, 40, 50, 60, 70, 80, Inf),
    #   labels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"),
    #   right = TRUE
    # ),
    # factor(age10, levels = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"))
    # ) %>% 
    mutate(imd5_origin = case_when(e2019_imd_10 ==1L | e2019_imd_10 ==2L ~ "1", 
                                   e2019_imd_10 ==3L | e2019_imd_10 ==4L ~ "2",
                                   e2019_imd_10 ==5L | e2019_imd_10 ==6L ~ "3",
                                   e2019_imd_10 ==7L | e2019_imd_10 ==8L ~ "4",
                                   e2019_imd_10 ==9L | e2019_imd_10 ==10L ~ "5"
    )) %>% 
    mutate(
      imd5_imp2 = {
        tbl <- table(imd5_imp, useNA = "no")
        if_else(is.na(imd5_imp), names(tbl)[which.max(tbl)], imd5_imp)
      }
    )
  
  ugi_symp6 <- ugi_symp5 %>% left_join(death)
  
  # Add end of follow-up date
  # could be due to death, or until 2021-03-29
  # see ONS death registration document
  ugi_symp6 <- ugi_symp6 %>% 
    mutate(
      death = if_else(!is.na(death_date), 1, 0),
      fu_end_date = if_else(!is.na(death_date), death_date, as.Date("2021-03-29")))
  
  return(ugi_symp6)
}

category_symp <- function(ugi_symp, symp_list, ng12_symp, red_flag){
  
  # Create indicator columns
  ugi_symp <- ugi_symp %>%
    mutate(
      # Check if symptoms_28day contains any NG12 symptom
      has_ng12_symp = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(ng12_symp, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains any red flag symptom
      has_red_flag = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(red_flag, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains any symptoms in the list
      has_any_symp = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(symp_list, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains weight loss
      has_weight_loss = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "weight_loss"))
      ),
      
      # Check if symptoms_28day contains weight_loss AND at least one other NG12 symptom
      has_weight_loss_plus = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        {
          # Create vector of NG12 symptoms excluding weight_loss
          other_ng12 <- ng12_symp[ng12_symp != "weight_loss"]
          # Check for weight_loss AND at least one other NG12 symptom
          as.integer(str_detect(symptoms_28day, "weight_loss") & 
                       str_detect(symptoms_28day, paste(other_ng12, collapse = "|")))
        }
      )
    ) %>% 
    
    # Create categorical variable to indicate risk
    mutate(
      
      risk_level = case_when(has_red_flag == 1 | (age_index >=55 & has_weight_loss_plus==1) ~ "has red flag symptoms",
                             has_ng12_symp == 1 ~ "has NG12 listed symptoms",
                             TRUE ~ "has extended symptoms only"
                             ),
      risk_level = factor(risk_level, levels = c("has extended symptoms only", 
                                                    "has NG12 listed symptoms",
                                                    "has red flag symptoms"
                                                    ))
    )
  
  return(ugi_symp)
}


category_symp_ng12 <- function(ugi_symp){
  
  # Code lists for "any" symptoms in NG12
  # Pancreatic cancer level 2 symptoms
  # combine with aged 60 + weight loss
  panc_symp_level2 <- c("diarrhoea", "painback", "epigastric_pain", "nausea/vomiting", 
                 "constipation", "new diabetes T1/2")
  
  # oesophageal and stomach level 2 symptoms
  # combine with age 55 + weight loss
  oeso_stom_symp_level21 <- c("epigastric_pain", "heartburn", "dyspepsia", "tr_dyspepsia")
  
  # combine with age 55 + raised platelet count
  oeso_stom_symp_level22 <- c("nausea/vomiting", "weight_loss", "heartburn", 
                              "dyspepsia", "tr_dyspepsia", "epigastric_pain")
  
  # combine with age 55 + nausea/vomitting
  oeso_stom_symp_level23 <- c("weight_loss", "heartburn", 
                              "dyspepsia", "tr_dyspepsia", "epigastric_pain")
  
  ugi_symp_ng12 <- ugi_symp %>% 
    mutate(
      
      # dysphagia
      has_dysphagia = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "dysphagia"))
      ),
      
      # haematemesis
      has_haematemesis = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "haematemesis"))
      ),
      
      # upper_abdo_mass
      has_upper_abdo_mass = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "upper_abdo_mass"))
      ),
      
      # jaundice
      has_jaundice = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "jaundice"))
      ),
      
      # weight loss
      has_weight_loss = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "weight_loss"))
      ),
      
      # treatment resistant dyspepsia
      has_tr_dyspepsia = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "tr_dyspepsia"))
      ),
      
      # hb_low
      has_hb_low = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "hb_low"))
      ),
      
      # plt_high
      has_plt_high = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "plt_high"))
      ),
      
      # epigastric_pain
      has_epigastric_pain = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "epigastric_pain"))
      ),
      
      # nausea/vomiting
      has_nausea = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, "nausea/vomiting"))
      ),
      
      # Check if symptoms_28day contains any symptoms in the list
      has_any_panc_symp_level2 = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(panc_symp_level2, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains any symptoms in the list
      has_any_oeso_stom_symp_level21 = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(oeso_stom_symp_level21, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains any symptoms in the list
      has_any_oeso_stom_symp_level22 = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(oeso_stom_symp_level22, collapse = "|")))
      ),
      
      # Check if symptoms_28day contains any symptoms in the list
      has_any_oeso_stom_symp_level23 = ifelse(
        is.na(symptoms_28day) | symptoms_28day == "",
        0,
        as.integer(str_detect(symptoms_28day, paste(oeso_stom_symp_level23, collapse = "|")))
      )
    ) %>% 
    mutate(
      
      ng12_red_flag = case_when(age_index >=40 & has_jaundice==1 ~ "panc 2ww",
                                
                                has_upper_abdo_mass == 1 ~ "stom/galb 2ww", 
                                
                                has_dysphagia == 1 ~ "oeso/stom 2ww",
                                
                                age_index >=55 & has_weight_loss==1 &
                                  has_any_oeso_stom_symp_level21 == 1 ~ "oeso/stom 2ww",
                                
                                age_index >=60 & has_weight_loss==1 & 
                                  has_any_panc_symp_level2 == 1 ~ "panc CT/USS",
                                
                                has_haematemesis == 1 ~ "oeso/stom endoscopy",
                                
                                age_index >=55 & has_tr_dyspepsia==1 ~ "oeso/stom endoscopy",
                                
                                age_index >=55 & has_hb_low==1 & 
                                  has_epigastric_pain == 1 ~ "oeso/stom endoscopy",
                                
                                age_index >=55 & has_plt_high==1 &
                                  has_any_oeso_stom_symp_level22 == 1 ~ "oeso/stom endoscopy",
                                
                                age_index >=55 & has_nausea==1 &
                                  has_any_oeso_stom_symp_level23 == 1 ~ "oeso/stom endoscopy",
                                
                                TRUE ~ "No red flag"
      )
    )
  
  return(ugi_symp_ng12)
}



















