# Create low haemoglobin binary variable
rm(list = ls())
library(haven)
library(tidyverse)
library(dtplyr)
library(data.table)
source("/data/WIPH-CanDetect/HealthEco/route.R")
source(file.path(scr, "data_clean", "fn_clean.R"))

# Load data ---------------------------------------------------------------
# Load unit code list for reference
# numunit <- read.delim(file.path(candetect, "documentation", "lookups", 
#                                 "aurum", "NumUnit.txt"))
# The original unit code list has a few errors after reading into R
# I used Excel to read it and convert into .csv file, which addressed this issue.
numunit <- read.csv(file.path(wd, "code_list", "NumUnit_correct.csv"))

# Load all haemoglobin medcodeids
hb_medcode <- read_xlsx(file.path(candetect, "codelists", 
                                 "CPRD_Aurum_Codelist_Library", "Lab_Tests",
                                 "cprd_aurum_hb.xlsx"))

# Load CPRD observations 2012-04-01 - 2025-01-01
# all medcodeid within "WIPH-CanDetect/codelists/CPRD_Aurum_Codelist_Library/Lab_tests/CPRD_Aurum_hb.dta")
cprd_hb <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_hb.rds"))

# Extract all used unit codes in order of frequency
# there are 84 unit codes here
# many of them are only used 1-10 times
unitcodes <- as.integer(names(sort(table(cprd_hb$numunitid), decreasing = T)))
# Look up the unit descriptions
# Export to manual check
unit_dic <- numunit[match(unitcodes, numunit$numunitid),]
unit_dic$freq <- sort(table(cprd_hb$numunitid), decreasing = T)
write.csv(unit_dic, file.path(output, "hb_unit.csv"))

# Load patient level data for sex identification
cprd_pat <- readRDS(file.path(wd, "create_study_pop_interm", "cprd_pat.rds")) %>% 
  mutate(e_patid = as.character(e_patid))

cprd_hb <- cprd_hb %>% 
  left_join(cprd_pat[, c("e_patid", "male")], by = "e_patid")

# Recode with plausible ranges --------------------------------------------

# 139 g/L has 43.7m observation
# 138 g/dL has 5.1m observation
# 17 units have 10-384 observations
# manually code them by checking with values against the plausible range
unit_dic2 <- unit_dic %>% 
  mutate(
    unit_cor = case_when(Description == "g/L" ~ "g/L",
                         Description == "g/dL" ~ "g/dL",
                         Description == "gm/dl" ~ "g/dL",
                         Description == "g/dL." ~ "g/dL",
                         Description == "%" ~ "unknow",
                         Description == "g/litre" ~ "g/L",
                         Description == "mmol/L" ~ "mol",
                         Description == "MG/DL" ~ "g/dL",
                         Description == "g/dlg/dl" ~ "g/dL",
                         Description == "gm/l" ~ "g/L",
                         Description == "mmol/mol" ~ "mol",
                         Description == "SI units" ~ "g/L",
                         Description == "l" ~ "g/L",
                         Description == "10*9/L" ~ "g/L",
                         Description == "g/d" ~ "g/dL",
                         Description == "l/l" ~ "g/L",
                         Description == "u/L" ~ "g/L",
                         Description == "g/mol" ~ "g/L",
                         Description == "ug/L" ~ "g/L",
                         grepl("mol", Description) ~ "mol",
                         TRUE ~ "unknow"
                         )
  )

# according to Blood_test_ranges from Ben and Kirsten
g_dl_range <- c(1, 22)
g_l_range <- c(10, 220)
# function to check the value's distance to the median of the range
# scale up g/dL by 10
dist_med <- function(value, g_dl_range, g_l_range){
  
  dis_g_dl <- abs(value - median(g_dl_range))
  dis_g_l <- abs(value - median(g_l_range))
  
  return(c("dis_g_dl" =  dis_g_dl *10, "dis_g_l" = dis_g_l))
}
# dist_med(21, g_dl_range, g_l_range)
# the threshold is 21
threshold <- 21

# Correct unit types against values through 3 steps
# Could skip some of them if think it is inappropriate 
# 1: For unit_cor = "unknown" or NA, check value against the plausible range and assign them
cprd_hb2 <- cprd_hb %>% 
  left_join(unit_dic2[, c("numunitid", "unit_cor")], by = "numunitid") %>% 
  lazy_dt() %>% 
  mutate(
    unit_cor = case_when(
      (is.na(unit_cor) | unit_cor == "unknow") & 
        value >= g_dl_range[1] & value < threshold ~ "g/dL",
      (is.na(unit_cor) | unit_cor == "unknow") & 
        value >= threshold & value <= g_l_range[2] ~ "g/L",
      TRUE ~ unit_cor
    )
  ) %>% 

# 2: For unit_cor = "mol", most values are out of the plausible range of mmol/L
# change them to g/dL or g/L against the plausible range
  mutate(
    unit_cor = case_when(
      unit_cor == "mol" & 
        value >= g_dl_range[1] & value < threshold ~ "g/dL",
      unit_cor == "mol" & 
        value >= threshold & value <= g_l_range[2] ~ "g/L",
      TRUE ~ unit_cor
    )
  ) %>% 

# 3: For unit_cor = "g/dL" or "g/L", check value against the plausible range
# keep them are, change type or code as impossible
  mutate(
    unit_cor = case_when(
      # keep those in the plausible range as they are
      unit_cor == "g/dL" & value >= g_dl_range[1] & value<= g_dl_range[2] ~ "g/dL",
      unit_cor == "g/L" & value >= g_l_range[1] & value<= g_l_range[2] ~ "g/L",
      # assign those to the other type
      # If the value is out of the plausible range of one, but inside the plausible range of the other
      unit_cor == "g/dL" & value > g_dl_range[2] & value<= g_l_range[2] ~ "g/L",
      unit_cor == "g/L" & value > g_dl_range[1] & value< g_l_range[1] ~ "g/dL",
      is.na(unit_cor) ~ NA,
      # code others "impossible"
      TRUE ~ "impossible"
    )
  ) %>% 

# convert g/L to g/dL
  mutate(
    value_gl = case_when(unit_cor == "g/L" ~ value,
                          unit_cor == "g/dL" ~ value*10,
                          TRUE ~ NA # where all values are out of both plausible ranges
                          )
  ) %>% 
  as.data.frame()

# Create the binary variable low_hb ---------------------------------------
hb_normal_m <- c(132, 166)
hb_normal_f <- c(116, 150)

cprd_hb2 <- cprd_hb2 %>% 
  mutate(
    hb_level = case_when(
      male = 1 & value_gl < hb_normal_m[1] ~ "low",
      male = 1 & value_gl >= hb_normal_m[1] & value_gl <=hb_normal_m[2] ~ "normal",
      male = 1 & value_gl > hb_normal_m[2] ~ "high",
      male = 0 & value_gl < hb_normal_f[1] ~ "low",
      male = 0 & value_gl >= hb_normal_f[1] & value_gl <=hb_normal_f[2] ~ "normal",
      male = 0 & value_gl > hb_normal_f[2] ~ "high",
      TRUE ~ "unknown"
    )
  )

# For "unknow", look up medcodeid
# refer to terms in hb_medcode from WS1
x <- table(cprd_hb2$medcodeid, useNA = "ifany")
medcode <- names(sort(x, decreasing = T))
rt <- hb_medcode[match(medcode, hb_medcode$medcodeid), c("medcodeid", "term")]
rt$freq <- sort(x, decreasing = T)
print(rt, n = nrow(rt))

# identify medcodeid linking to low haemoglobin
low_medcode <- rt[grepl("low", rt$term), c("medcodeid")] # below is included
high_medcode <- rt[grepl("high|Increased|above", rt$term), c("medcodeid")]
normal_medcode <- rt[grepl("\\bnormal\\b|within", rt$term), c("medcodeid")]
abnormal_medcode <- rt[grepl("Abnormal|outside", rt$term), c("medcodeid")]

cprd_hb2 <- cprd_hb2 %>% 
  # mutate(medcodeid = as.character.integer64(medcodeid)) %>% 
  mutate(
    hb_level2 = case_when(
      hb_level == "unknown" & 
        medcodeid %in% as.character(low_medcode$medcodeid) ~ "low",
      hb_level == "unknown" & 
        medcodeid %in% as.character(high_medcode$medcodeid) ~ "high",
      hb_level == "unknown" & 
        medcodeid %in% as.character(normal_medcode$medcodeid) ~ "normal",
      hb_level == "unknown" & 
        medcodeid %in% as.character(abnormal_medcode$medcodeid) ~ "abnormal",
      TRUE ~ hb_level
    )
  ) %>% 
  #         abnormal     high      low   normal  unknown
  # high           0   819475        0        0        0
  # low            0        0 20343490        0        0
  # normal         0        0        0 28042610        0
  # unknown      880     3500    26250     1842   254296
  select(e_patid, obsid, date, hb_level2) %>% 
  rename(haemoglobin = hb_level2)

saveRDS(cprd_hb2, file.path(wd, "create_study_pop_interm", "haemoglobin_clean.rds")) 
# consider assign abnormal to low

