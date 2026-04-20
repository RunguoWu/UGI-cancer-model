# Create death indicator
rm(list = ls())

library(tidyverse)
library(data.table)
library(dtplyr)
library(bit64)
library(skimr)

source("/data/WIPH-CanDetect/HealthEco/route.R")

# linked death records
death <- fread(file.path(link_raw, "e_aurum_death_patient_23_002840_dm.txt"))

# code if a cancer type is the cause of death
death$death_oeso <- 0
death$death_panc <- 0
death$death_galb <- 0
death$death_stom <- 0
death$death_upGI <- 0
death$death_cancer <- 0

for (i in c("", 1:15)) {
  text <- paste0("cause", i) 
  death$death_oeso[grepl("^C15", death[[text]])] <- 1 
  death$death_panc[grepl("^C25", death[[text]])] <- 1
  death$death_stom[grepl("^C16", death[[text]])] <- 1
  death$death_galb[grepl("^C23|^C24", death[[text]])] <- 1
  death$death_upGI[grepl("^C15|^C16|^C23|^C24|^C25", death[[text]])] <- 1
  death$death_cancer[grepl("^C", death[[text]]) & !grepl("^C44", death[[text]])] <- 1
}

death <- death %>% 
  lazy_dt() %>% 
  mutate(dor = as.Date(dor, "%d/%m/%Y"), 
         dod = as.Date(dod, "%d/%m/%Y")) %>% 
  mutate(death_date = if_else(is.na(dod), dor, dod)) %>% 
  mutate(e_patid=as.character.integer64(e_patid)) %>% 
  select(e_patid, death_date, 
         death_oeso, death_panc, death_stom, 
         death_galb, death_upGI, death_cancer) %>% 
  as.data.frame()

saveRDS(death, file.path(wd, "create_study_pop_interm", "death.rds"))
