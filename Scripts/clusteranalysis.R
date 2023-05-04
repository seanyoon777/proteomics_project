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
load_lib("ggrepel")

get_biodata <- function(path) {
  dir <- "/Users/seonghyunyoon/Downloads/proteomics_project/data"
  read.csv(paste(dir, path, sep = "/"))
}

plasma_meta <- get_biodata("Proteomic/Plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
plasma_prot <- get_biodata("Proteomic/Plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
plasma_patient <- get_biodata("Proteomic/Plasma/Plasma_BasicDemo.csv")
CSF_meta <- get_biodata("Proteomic/CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
CSF_prot <- get_biodata("Proteomic/CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
CSF_patient <- get_biodata("Proteomic/CSF/CSF_BasicDemo.csv")

patient_meta <- get_biodata("Proteomic/commonfile/WashU_ADNI_demographic.csv")

all_prots_ID <- names(plasma_prot)[4:ncol(plasma_prot)] 

protMeta <- read_csv("/Users/seonghyunyoon/Downloads/proteomics_project/data/ProteinMetadata_with_LOD.csv", col_types = )
protMeta$SeqId = str_replace(protMeta$SeqId, "-", ".")
protMeta <- filter(protMeta, (Organism == "Human" | Organism == "HIV-1" | Organism == "HIV-2") & Type == "Protein")
all_prots <- str_remove(all_prots_ID, "^X")
my_order <- match(protMeta$SeqId, all_prots)
all_prots <- protMeta$Key_2
nprots <- length(all_prots)

plasma_prot <- plasma_prot[c(1:3, 3 + my_order)]
colnames(plasma_prot) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)

CSF_prot <- CSF_prot[c(1:3, 3 + my_order)]
colnames(CSF_prot) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)


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
  
# strict: set date window
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



# lenient time window 
date_window <- 2000
patient_lenient <- patientdata_summary[patientdata_summary$drawdate_diff <= date_window, ]
intake_lenient <- intake_fil[patientdata_summary$drawdate_diff <= date_window, ]

intake_lenient_lm_models <- list()
for(i in 1:length(all_prots)) {
  intake_lenient_lm_models[[i]] <- summary(lm(intake_lenient[, i] ~ 
                                              relevel(factor(patient_lenient$final_status), ref = "CO") + 
                                              patient_lenient$Sex +
                                              patient_lenient$age + 
                                              patient_lenient$storage_days))
}



intake_lenient_lm_summary <- data.frame()
for (i in 1:length(all_prots)) {
  intake_lenient_tidy_result <- data.frame(intake_lenient_lm_models[[i]]$coefficients) %>% 
    select(Estimate, Std..Error, Pr...t..)
  colnames(intake_lenient_tidy_result) <- c("Coefficient", "StdError", "Pval")
  intake_lenient_tidy_result <- intake_lenient_tidy_result[2:nrow(intake_lenient_tidy_result), ] %>% 
    mutate(xVar = c("AD", "FTD", "Non-AD dementia", "Male", "Age", "Storage_days"), 
           Protein = all_prots[i]) %>% 
    select(Protein, xVar, Coefficient, StdError, Pval)
  rownames(intake_lenient_tidy_result) <- c()
  
  intake_lenient_lm_summary <- rbind(intake_lenient_lm_summary, intake_lenient_tidy_result)
}



generate_lmsummary <- function(data, date_window, data_type = c("ADRC", "Knight"))


generate_volcanodata <- function(data, variable) {
  volcano_data <- data %>% filter(xVar == variable) %>% 
  #group_by(Protein) %>%
  mutate(padj = p.adjust(Pval, method = "fdr")) %>%
  mutate(qval = -log10(padj), 
         log2fc = Coefficient) %>% 
  select(Protein, Pval, padj, qval, log2fc) %>% 
  mutate(diffexpressed = ifelse(log2fc > 0 & qval >= 1,
                                yes = "Upregulated", 
                                no = ifelse(log2fc < 0 & qval >= 1,
                                            yes = "Downregulated", no="none")))
  return(volcano_data)
}

generate_volcanoplot <- function(volcano_data, variable) {
  ggplot(volcano_data, aes(x = log2fc, y = qval, color = factor(diffexpressed))) + 
    geom_point(size = 0.5, alpha = 0.8, na.rm = T) +
    #geom_text_repel(max.overlaps = 10, aes(label = delabel)) + 
    theme_bw(base_size = 16) +
    theme(legend.position = "none") +
    #ggtitle(label = paste(str_split(plot_title, '_', simplify = T), sep = "", collapse = " ")) + 
    xlab("log2 FC") +
    ylab(expression(-log[10]("q"))) +
    scale_color_manual(values = c("Upregulated" = "indianred1",
                                  "Downregulated" = "royalblue1",
                                  "none" = "grey60")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
          plot.background = element_rect(fill = "white"),
          axis.title = element_text(size = 9.5),
          axis.text = element_text(size = 7),
          plot.title = element_text(size = 10, hjust = 0.5), 
          axis.line = element_blank()) + 
    annotate("text", x = min(volcano_data$log2fc)*1.1, y = max(volcano_data$qval), 
             label = variable, hjust = 0, size = 4)
}


data_knight <- c(intake_strict_lm_summary, intake_lenient_lm_summary)
xVars <- c("AD", "FTD", "Non-AD dementia", "Male", "Age", "Storage_days")

knight_strict_volcanodata <- list()
for (i in 1:length(xVars)) {
  knight_strict_volcanodata[[i]] <- generate_volcanodata(intake_strict_lm_summary, 
                                                 variable = xVars[i])
}
knight_strict_volcanoplot <- list()
for (i in 1:length(xVars)) {
  knight_strict_volcanoplot[[i]] <- generate_volcanoplot(knight_strict_volcanodata[[i]], 
                                                         variable = xVars[i])
}

gridExtra::grid.arrange(grobs = knight_strict_volcanoplot, ncol = 3)


# 1. confirm if group_by is right? 
# 2. generate new volcano plots 
# 3. Send the new volcanodata datasets 
# 4. Clustering 







#This time for lenient datewindow
knight_lenient_volcanodata <- list()
for (i in 1:length(xVars)) {
  knight_lenient_volcanodata[[i]] <- generate_volcanodata(intake_lenient_lm_summary, 
                                                         variable = xVars[i])
}
knight_lenient_volcanoplot <- list()
for (i in 1:length(xVars)) {
  knight_lenient_volcanoplot[[i]] <- generate_volcanoplot(knight_lenient_volcanodata[[i]], 
                                                         variable = xVars[i])
}

gridExtra::grid.arrange(grobs = knight_lenient_volcanoplot, ncol = 3)


# export some data
(knight_strict_volcanodata[[1]])$Protein <- all_prots


# same thing with stanford data?
barcodes2 <- get_biodata("CSF_plasma_matched_barcodes.csv")
CSF_patient2 <- get_biodata("CSF_metadata_2023-03-04.csv")
CSF_protein2 <- get_biodata("CSFProts.log10.noLODFilter.csv")
plasma_patient2 <- get_biodata("Plasma_metadata_samplesOnly_2023-03-04.csv")
plasma_protein2 <- get_biodata("plasmaProts.log10.csv")

CSF_patient2$Age <- gsub('^"&"$', '', CSF_patient2$Age)
CSF_patient2$Age <- as.numeric(CSF_patient2$Age)

common_prots2 <- intersect(names(CSF_protein2), names(plasma_protein2))
CSF_protein2 <- data.frame(Barcode = CSF_patient2$Barcode, CSF_protein2[, common_prots2]) 
plasma_protein2 <- data.frame(Barcode = plasma_protein2$Barcode, plasma_protein2[, common_prots2])


CSF_protein2_fil <- CSF_protein2[CSF_protein2$Barcode %in% barcodes2$CSF_Barcode, ]
CSF_patient2_fil <- CSF_patient2[CSF_patient2$Barcode %in% barcodes2$CSF_Barcode, ]
plasma_protein2_fil <- plasma_protein2 %>%
  filter(Barcode %in% barcodes2$Plasma_Barcode) %>% 
  arrange(match(Barcode, barcodes2$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, barcodes2$Plasma_Barcode)))

age2_filtered <- plasma_patient2 %>% 
  select(Age, Barcode) %>% 
  rename(Plasma_Barcode = Barcode) %>% 
  filter(Plasma_Barcode %in% barcodes2$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes2$Plasma_Barcode)) %>% 
  mutate(CSF_Barcode = CSF_patient2_fil$Barcode) %>% 
  select(Age, CSF_Barcode, Plasma_Barcode)

patientdata_filtered <- plasma %>% 
  rename(Plasma_Barcode = Barcode, Plasma_Zscore = ConnectivityZscore) %>% 
  filter(Plasma_Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes$Plasma_Barcode)) %>% 
  mutate(CSF_Zscore = CSF_patientdata_filtered$ConnectivityZscore, 
         CSF_Barcode = CSF_patientdata_filtered$Barcode) %>% 
  .[, c(ncol(.), 1:(ncol(.) - 1))] %>% 
  .[, c(1:(ncol(.) - 2), ncol(.), ncol(.) - 1)]





# CLUSTERING ANALYSIS
intake_fil












