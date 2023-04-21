rm(list = ls())

#setting foundations
load_lib <- function(packages, repos = "http://cran.us.r-project.org") {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = repos)
    }
    library(package, character.only = TRUE)
  }
}

load_lib("tidyverse")
load_lib("caret")
load_lib("data.table")
load_lib("corrr")
load_lib("igraph")
load_lib("cluster")
load_lib("factoextra")
load_lib("purrr")
load_lib("corrplot")
load_lib("dplyr")
load_lib("ggplot2")
load_lib("pls")
load_lib("ggraph")
load_lib("circlize")
load_lib("Cairo")
load_lib("ComplexHeatmap")
load_lib("loessclust")
load_lib("broom")

get_biodata <- function(path) {
  dir <- "/Users/seonghyunyoon/Downloads/proteomics_project/data/Proteomic"
  read.csv(paste(dir, path, sep = "/"))
}

plasma_meta <- get_biodata("Plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
plasma_prot <- get_biodata("Plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
plasma_patient <- get_biodata("Plasma/Plasma_BasicDemo.csv")
CSF_meta <- get_biodata("CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
CSF_prot <- get_biodata("CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
CSF_patient <- get_biodata("CSF/CSF_BasicDemo.csv")

patient_meta <- get_biodata("commonfile/WashU_ADNI_demographic.csv")

nprots <- ncol(plasma_prot)
all_prots <- names(plasma_prot)[4:nprots]

# Logarithmic
plasma_prot[4:nprots] <- log10(plasma_prot[4:nprots])
CSF_prot[4:nprots] <- log10(CSF_prot[4:nprots])

## 0.1. Wrangle CSF and plasma dataset 
match_by <- function(data, refdata) {
  data %>% 
    filter(PA_DB_UID %in% refdata$PA_DB_UID) %>% 
    arrange(match(PA_DB_UID, refdata$PA_DB_UID)) %>% 
    slice(order(match(PA_DB_UID, patient_meta$PA_DB_UID)))
}

plasma_patient <- plasma_patient %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed AD", 
                   "AD", final_cc_status.updated)) %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed Control", 
                   "CO", final_cc_status.updated)) %>% 
  filter(final_cc_status.updated %in% c("AD", "CO")) %>% 
  mutate(drawdate = strptime(drawdate, format = "%m/%d/%y"), 
         DateOfBirth = gsub("^\\d{2}", "19\\", strptime(DateOfBirth, format = "%m/%d/%y")))

plasma_meta <- match_by(plasma_meta, plasma_patient)
plasma_prot <- match_by(plasma_prot, plasma_patient)
plasma_prot_zscored <- plasma_prot
plasma_prot_zscored[all_prots] <- scale(plasma_prot[all_prots])

CSF_patient <- CSF_patient %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed AD", 
                   "AD", Last.status)) %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed Control", 
                   "CO", Last.status)) %>% 
  filter(Last.status %in% c("AD", "CO")) %>% 
  mutate(DrawDate = strptime(DrawDate, format = "%m/%d/%y"))

## 0.2. Wrangle intake dataset ----
common_ID <- intersect(CSF_patient$PA_DB_UID, plasma_patient$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .)

CSF_meta <- match_by(CSF_meta, CSF_patient)
CSF_prot <- match_by(CSF_prot, CSF_patient)
CSF_prot_zscored <- CSF_prot
CSF_prot_zscored[all_prots] <- scale(CSF_prot[all_prots])
CSF_meta_fil <- match_by(CSF_meta, common_ID)
CSF_prot_fil <- match_by(CSF_prot, common_ID)
CSF_patient_fil <- match_by(CSF_patient, common_ID)
plasma_meta_fil <- match_by(plasma_meta, common_ID)
plasma_prot_fil <- match_by(plasma_prot, common_ID)
plasma_patient_fil <- match_by(plasma_patient, common_ID)
patient_meta_fil <- match_by(patient_meta, common_ID)

intake_fil <- 1 - CSF_prot_fil[4:nprots] / plasma_prot_fil[4:nprots]
patientdata_fil <- merge(CSF_patient_fil, plasma_patient_fil, by = "PA_DB_UID") %>% 
  filter(Last.status %in% c("AD", "CO")) %>% 
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status) %>%
  mutate(drawdate_diff = abs(difftime(Plasma_drawdate, CSF_drawdate, units = "days")))

intake_fil <- intake_fil[patientdata_fil$Plasma_drawstatus == 
                           patientdata_fil$CSF_drawstatus,  ]
patientdata_fil <- patientdata_fil[patientdata_fil$Plasma_drawstatus == 
                                     patientdata_fil$CSF_drawstatus, ] %>% 
  mutate(drawage = (Plasma_drawage + CSF_drawage) / 2)

## 0.2. Wrangle common patientdata
drawdate_diff <- difftime(CSF_patient_fil$DrawDate, plasma_patient_fil$drawdate, units = "days") %>%
  data.frame(Diff = ., absDiff = abs(.))







