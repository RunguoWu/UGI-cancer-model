# Clean CPRD patient level data #

rm(list = ls())
library(parallel)
library(doSNOW)
library(tidyverse)
source("/data/WIPH-CanDetect/HealthEco/route.R")

## Combine CPRD patient data
cl <- makePSOCKcluster(4)
clusterExport(cl, c("cprd_raw"))
clusterEvalQ(cl, library(tidyverse))
registerDoSNOW(cl)

cprd_pat <- foreach(i = 0:43) %dopar% { # patients are allocated to 44 files
  
  file_name <- paste0("e_aurum_patlist", i, "_extract_patient_001.txt")
  
  x <- read.delim(file.path(cprd_raw, "patient_txt", file_name), header = T)
  
  x <- x %>% 
    mutate(dob_imp = as.Date(paste0(yob, "-06-30")),
           regstartdate = as.Date(regstartdate, format = "%d/%m/%Y"),
           regenddate = as.Date(regenddate, format = "%d/%m/%Y"),
           cprd_ddate = as.Date(cprd_ddate, format = "%d/%m/%Y"),
           male = ifelse(gender==1, 1, ifelse(gender==2, 0, NA))
    ) %>% 
    select(e_patid, e_pracid, gender, male, dob_imp, cprd_ddate,
           regstartdate, regenddate, acceptable)
  
  return(x)         
}
stopCluster(cl)

cprd_pat <- do.call("rbind", cprd_pat)
cprd_pat$acceptable <- NULL # all acceptable = 1

cprd_pat$e_patid <- as.character(cprd_pat$e_patid)

saveRDS(cprd_pat, file.path(wd, "create_study_pop_interm", "cprd_pat.rds"))


# Linked data -------------------------------------------------------------
## IMD----
file_name <- paste0("e_aurum_patient_2019_imd_23_002840.txt")
imd <- read.delim(file.path(link_raw, file_name), header = T)

# very few IMD data is missing
# Impute using the most frequent IMD10 value within the same GP practice 
imd2 <- imd %>% 
  lazy_dt() %>% 
  group_by(e_pracid) %>% 
  # Calculate mode excluding NAs
  mutate(
    practice_mode = {
      tbl <- table(e2019_imd_10, useNA = "no")  # Exclude NAs
      if(length(tbl) > 0) {
        as.integer(names(tbl)[which.max(tbl)])
      } else {
        NA_real_
      }
    },
    imd10_imp = if_else(is.na(e2019_imd_10), practice_mode, e2019_imd_10)
  ) %>% 
  ungroup() %>% 
  select(-practice_mode) %>% 
  as.data.frame()

imd2 <- imd2 %>% 
  # 1 = least deprived
  mutate(imd5_imp = case_when(imd10_imp ==1L | imd10_imp ==2L ~ "1", 
                          imd10_imp ==3L | imd10_imp ==4L ~ "2",
                          imd10_imp ==5L | imd10_imp ==6L ~ "3",
                          imd10_imp ==7L | imd10_imp ==8L ~ "4",
                          imd10_imp ==9L | imd10_imp ==10L ~ "5"
  ))

imd2$e_patid <- as.character(imd2$e_patid)

saveRDS(imd2, file.path(wd, "create_study_pop_interm", "imd.rds"))


## Read in HES patient data ----
hes_pat <- read.delim(file.path(dt_link, "e_aurum_hes_patient_23_002840_dm.txt"),
                      header = T)