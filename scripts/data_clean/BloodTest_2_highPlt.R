# Create raised platelet binary variable
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
plt_medcode <- read_xlsx(file.path(candetect, "codelists", 
                                 "CPRD_Aurum_Codelist_Library", "Lab_Tests",
                                 "cprd_aurum_platelets.xlsx"))

# Load CPRD observations 2012-04-01 - 2025-01-01
# all medcodeid within "WIPH-CanDetect/codelists/CPRD_Aurum_Codelist_Library/Lab_tests/CPRD_Aurum_hb.dta")
cprd_plt <- readRDS(file = file.path(wd, "create_study_pop_interm", "cprd_obs_plt.rds"))

# Extract all used unit codes in order of frequency
# there are 84 unit codes here
# many of them are only used 1-10 times
unitcodes <- as.integer(names(sort(table(cprd_plt$numunitid), decreasing = T)))
# Look up the unit descriptions
# Export to manual check only
unit_dic <- numunit[match(unitcodes, numunit$numunitid),]
unit_dic$freq <- sort(table(cprd_plt$numunitid), decreasing = T)
# write.csv(unit_dic, file.path(output, "plt_unit.csv"))


# Create high plt using the range -----------------------------------------

# Consider all unit 10^9/L, regardless the unit label, as we don't have an alternative
# Check them against the plausible range
plt_range <- c(10, 1000)
plt_normal <- c(150, 400)

cprd_plt2 <- cprd_plt %>% 
  mutate(
    plt_level = case_when(
      value >= plt_range[1] & value < plt_normal[1] ~ "low",
      value >=plt_normal[1] & value <=plt_normal[2] ~ "normal",
      value > plt_normal[2] & value <= plt_range[2] ~ "high",
      TRUE ~ "unknown")
  )


# For "unknow", look up medcodeid
# refer to terms in hb_medcode from WS1
x <- table(cprd_plt2$medcodeid, useNA = "ifany")
medcode <- names(sort(x, decreasing = T))
rt <- plt_medcode[match(medcode, plt_medcode$medcodeid), c("medcodeid", "term")]
rt$freq <- sort(x, decreasing = T)
print(rt, n = nrow(rt))

# identify medcodeid linking to low, high, normal or abnormal plt
# no medcodeid for low plt
high_medcode <- rt[grepl("Increased", rt$term), c("medcodeid")]
normal_medcode <- rt[grepl("count normal|within", rt$term), c("medcodeid")]
abnormal_medcode <- rt[grepl("count abnormal|outside", rt$term), c("medcodeid")]

cprd_plt2 <- cprd_plt2 %>% 
  mutate(
    plt_level2 = case_when(
      plt_level == "unknown" & 
        medcodeid %in% as.character(high_medcode$medcodeid) ~ "high",
      plt_level == "unknown" & 
        medcodeid %in% as.character(normal_medcode$medcodeid) ~ "normal",
      plt_level == "unknown" & 
        medcodeid %in% as.character(abnormal_medcode$medcodeid) ~ "abnormal",
      TRUE ~ plt_level
    )
  ) %>% 
  # tab plt_level plt_level2
  #         abnormal     high      low   normal  unknown
  # high           0  2373847        0        0        0
  # low            0        0  2234443        0        0
  # normal         0        0        0 44044249        0
  # unknown    16667        3        0     2343   282413
  select(e_patid, obsid, date, plt_level2) %>% 
  rename(platelet = plt_level2)

saveRDS(cprd_plt2, file.path(wd, "create_study_pop_interm", "platelet_clean.rds"))

