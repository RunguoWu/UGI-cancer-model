library(data.table)
library(bit64)

dt <- readRDS(file.path(wd, "study_pop_symp_stageImputed.rds"))
ncras <- fread(file.path(raw_dat, "aurum_ncras", "e_aurum_ncras_tumour_23_002840.txt"))

morph <- ncras %>% 
  filter(e_patid %in% dt$e_patid) %>% 
  mutate(cancerdate = as.Date(diagnosisdatebest, format = "%d/%m/%Y"), 
         e_patid = as.character.integer64(e_patid)
         ) %>% 
  filter(cancerdate >= "2014-04-01") %>% 
  mutate(icd10_o2_bhv = paste0(morph_icd10_o2, "/", behaviour_icd10_o2)) %>% 
  select(e_patid, cancerdate, icd10_o2_bhv, morph_icd10_o2, behaviour_icd10_o2) %>% 
  distinct()

dt2 <- dt %>% left_join(morph, by = c("e_patid", "cancerdate")) %>% 
  select(e_patid, cancerdate, site, icd10_o2_bhv)

icd10_o2_all <- unique(dt2$icd10_o2_bhv)


dt2$histology_group <- NA
# oeso ----
names(table(dt2[dt2$site=="oeso",]$icd10_o2_bhv,useNA = "ifany"))

scc_codes <- c("8051/3", "8070/2", "8070/3", "8070/5", "8071/3", "8072/3", "8074/3", "8090/3")

ac_codes <- c("8140/2", "8140/3", "8144/3", "8145/3", "8200/3", "8211/3", "8240/3",
              "8244/3", "8246/3", "8260/3", "8263/2", "8310/3", "8312/3",
              "8480/3", "8481/3", "8490/3", "8560/3")

# Assign category
dt2$histology_group[dt2$site=="oeso"] <-ifelse(dt2$icd10_o2_bhv[dt2$site=="oeso"] %in% scc_codes, "Squamous cell carcinoma",
                                               ifelse(dt2$icd10_o2_bhv[dt2$site=="oeso"] %in% ac_codes, "Adenocarcinoma",
                                                      "Other / Unclassified"))

# Check results
table(dt2$histology_group[dt2$site=="oeso"], useNA = "ifany")

# panc ----
names(table(dt2[dt2$site=="panc",]$icd10_o2_bhv,useNA = "ifany"))

# Define morphology groups
pdac_codes <- c("8010/3", "8012/3", "8020/3", "8021/3", "8022/3", "8033/3",
                "8140/3", "8141/3", "8145/3", "8210/3", "8430/3", "8452/3",
                "8470/3", "8471/3", "8480/3", "8481/3", "8490/3", "8500/3",
                "8550/3", "8560/3")

neuro_codes <- c("8151/3", "8154/3", "8240/3", "8244/3", "8245/3", "8246/3",
                 "8310/3", "8312/3")

# Assign histology group
dt2$histology_group[dt2$site=="panc"] <- ifelse(dt2$icd10_o2_bhv[dt2$site=="panc"] %in% pdac_codes, "Pancreatic ductal adenocarcinoma",
                                                ifelse(dt2$icd10_o2_bhv[dt2$site=="panc"] %in% neuro_codes, "Neuroendocrine pancreatic tumour",
                                                       "Other / Unclassified"))

# Check counts
table(dt2$histology_group[dt2$site=="panc"], useNA = "ifany")

# stom ----
names(table(dt2[dt2$site=="stom",]$icd10_o2_bhv,useNA = "ifany"))

adenoca_codes <- c("8140/2", "8140/3", "8142/3", "8144/3", "8145/3", "8160/3",
                   "8210/2", "8210/3", "8211/3", "8260/3", "8263/3", "8310/3", "8323/3",
                   "8480/3", "8481/3", "8490/3", "8560/3")

neuro_codes <- c("8240/3", "8240/6", "8244/3", "8246/3")

squamous_codes <- c("8070/2", "8070/3", "8071/3")

other_codes <- c("8000/3","8010/3","8012/3","8020/3","8032/3","8041/3","8045/3",
                 "8800/3","8890/3","8990/1","8990/3","9051/3","9680/3")

dt2$histology_group[dt2$site=="stom"] <- ifelse(dt2$icd10_o2_bhv[dt2$site=="stom"] %in% adenoca_codes, "Adenocarcinoma",
                                                ifelse(dt2$icd10_o2_bhv[dt2$site=="stom"] %in% neuro_codes, "Neuroendocrine tumour",
                                                       ifelse(dt2$icd10_o2_bhv[dt2$site=="stom"] %in% squamous_codes, "Squamous / Adenosquamous carcinoma",
                                                              ifelse(dt2$icd10_o2_bhv[dt2$site=="stom"] %in% other_codes, "Other / Unclassified", NA))))

table(dt2$histology_group[dt2$site=="stom"], useNA = "ifany")

# galb ----
names(table(dt2[dt2$site=="galb",]$icd10_o2_bhv,useNA = "ifany"))

# Define histology groups
adenoca_codes <- c("8140/2", "8140/3", "8144/3", "8160/3", "8162/3",
                   "8260/3", "8261/3", "8263/3", "8323/3",
                   "8480/3", "8481/3", "8490/3", "8500/3", "8503/3", "8560/3")

squamous_codes <- c("8050/3","8070/3")

neuro_codes <- c("8240/3","8244/3","8246/3")

other_codes <- c("8000/3","8010/3","8012/3","8020/3","8021/3","8041/3")

# Assign category
dt2$histology_group[dt2$site=="galb"] <- ifelse(dt2$icd10_o2_bhv[dt2$site=="galb"] %in% adenoca_codes, "Adenocarcinoma",
                                                ifelse(dt2$icd10_o2_bhv[dt2$site=="galb"] %in% squamous_codes, "Squamous / Adenosquamous carcinoma",
                                                       ifelse(dt2$icd10_o2_bhv[dt2$site=="galb"] %in% neuro_codes, "Neuroendocrine tumour",
                                                              ifelse(dt2$icd10_o2_bhv[dt2$site=="galb"] %in% other_codes, "Other / Unclassified", NA))))

table(dt2$histology_group[dt2$site=="galb"], useNA = "ifany")
