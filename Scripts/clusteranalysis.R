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
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status) %>%
  group_by(PA_DB_UID) %>% 
  mutate(drawdate_diff = abs(difftime(Plasma_drawdate, CSF_drawdate, units = "days")), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>% 
  mutate(storage_days = as.numeric(difftime(today(), avg_drawdate, units = "days"))) %>% 
  as.data.frame() 

patientdata_summary <- patientdata_fil %>% 
  mutate(final_status = ifelse(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus), 
         age = abs(Plasma_drawage + CSF_drawage) / 2) %>% 
  select(PA_DB_UID, Sex, age, drawdate_diff, avg_drawdate, final_status, storage_days) 

intake_fil <- intake_fil[patientdata_summary$final_status %in% c("AD", "CO", "FTD", "Non-AD dementia"), ]
patientdata_summary <- patientdata_summary %>% 
  filter(patientdata_summary$final_status %in% c("AD", "CO", "FTD", "Non-AD dementia"))
  
# set date window
date_window <- 180
patient_strict <- patientdata_summary[patientdata_summary$drawdate_diff <= date_window, ]
intake_strict <- intake_fil[patientdata_summary$drawdate_diff <= date_window, ]


for(i in 1:length(all_prots)) {
  intake_strict_lm_models[[i]] <- summary(lm(intake_strict[, i] ~ 
                                                relevel(factor(patient_strict$final_status), ref = "CO") + 
                                                patient_strict$Sex +
                                                patient_strict$age + 
                                                patient_strict$storage_days))
}


intake_strict_lm_summary <- data.frame()
for (i in 1:length(all_prots)) {
  intake_strict_tidy_result <- data.frame(intake_strict_lm_models[[i]]$coefficients) %>% 
    select(Estimate, Std..Error, Pr...t..)
  colnames(intake_strict_tidy_result) <- c("Coefficient", "StdError", "Pval")
  intake_strict_tidy_result <- intake_strict_tidy_result[2:nrow(intake_strict_tidy_result), ] %>% 
    mutate(xVar = c("AD", "FTD", "Non-AD dementia", "Male", "Age", "Storage_days"), 
           Protein = all_prots[i]) %>% 
    select(Protein, xVar, Coefficient, StdError, Pval)
  rownames(intake_strict_tidy_result) <- c()
  
  intake_strict_lm_summary <- rbind(intake_strict_lm_summary, intake_strict_tidy_result)
}

volcano_strict_AD <- intake_strict_lm_summary[intake_strict_lm_summary$xVar == "AD", ] %>% 
  mutate(qval = -log10(p.adjust(Pval, method = "fdr")), 
         log2fc = Coefficient) %>% 
  select(Protein, qval, log2fc)

ggplot(volcano_strict_AD, aes(x = log2fc, y = qval)) + 
  geom_point()



DEPlot <- ggplot(volcanoPlot, aes(x = log2fc, y = qvalue, color = factor(diffexpressed))) +
  geom_point(size = 3, alpha = 0.8, na.rm = T) + # add gene points
  geom_text_repel(max.overlaps = 10, aes(label = delabel)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  ggtitle(label = paste(str_split(plot_title, '_', simplify = T), sep = "", collapse = " ")) +
  xlab("log2 FC") +
  ylab(expression(-log[10]("q"))) +
  scale_color_manual(values = c("Upregulated" = "indianred1",
                                "Downregulated" = "royalblue1",
                                "none" = "grey60"))



## 0.2. Wrangle common patientdata
drawdate_diff <- difftime(CSF_patient_fil$DrawDate, plasma_patient_fil$drawdate, units = "days") %>%
  data.frame(Diff = ., absDiff = abs(.))






# set date window
date_window <- 180
patient_strict <- patientdata_fil[patientdata_fil$drawdate_diff <= date_window, ]
intake_strict <- intake_fil[patientdata_fil$drawdate_diff <= date_window, ]

# split dataset
intake_strict_AD <- intake_strict[patient_strict$Plasma_drawstatus == "AD", ]
intake_strict_CO <- intake_strict[patient_strict$Plasma_drawstatus == "CO", ]
split_intake_strict_age <- split(patient_strict$drawage, patient_strict$Plasma_drawstatus)
intake_strict_age_AD <- split_intake_strict_age$AD
intake_strict_age_CO <- split_intake_strict_age$CO

intake_strict_agemax <- min(max(intake_strict_age_AD), max(intake_strict_age_CO))
intake_strict_agemin <- max(min(intake_strict_age_AD), min(intake_strict_age_CO))
age_resolution <- 0.25
intake_age_seq <- seq(intake_strict_agemin, intake_strict_agemax, by = age_resolution)

loess_predict <- function(x, Y) {
  models <- vector("list", length(all_prots))
  for (i in 1:length(all_prots)) {
    models[[i]] <- loess(Y[, i] ~ x)
  }
  data <- data.frame(age = intake_age_seq)
  for (i in 1:length(all_prots)) {
    data[, i + 1] <- predict(models[[i]], intake_age_seq)
  }
  names(data) <- c("Age", all_prots)
  return(data)
}

models <- vector("list", length(all_prots))
for (i in 1:length(all_prots)) {
  models[[i]] <- loess(intake_strict_AD[, i] ~ intake_strict_age_AD)
}
data <- data.frame(age = intake_age_seq)
for (i in 1:length(all_prots)) {
  data[, i + 1] <- predict(models[[i]], intake_age_seq)
}

names(data) <- c("Age", all_prots)
intake_strict_AD_pred <- loess_predict(intake_strict_age_AD, intake_strict_AD)



# linear model, but in two ways: one with stringent dates, and one with up to 2000
# run lm on Knight, then create volcano plots for both stanford and the knight linear models 
intake_strict_lm_summary <- vector("list", length(all_prots))
for(i in 1:length(all_prots)) {
  intake_strict_lm_summary[[i]] <- summary(lm(intake_filtered[, i + 3] ~ 
                                       patient_strict$Plasma_drawstatus + 
                                       patient_strict$Sex +
                                       patient_strict$drawage))
}





