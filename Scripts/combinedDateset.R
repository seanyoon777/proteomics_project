source("Scripts/Helpers/init.R")
source("Scripts/Helpers/differentialAnalysis.R")
source("Scripts/Helpers/hierClust.R")
source("Scripts/Helpers/demoHist.R")

dir <- "/labs/twc/jarod/Data"
get_biodata <- function(path) {
  read.csv(paste(dir, path, sep = "/"))
}


# same thing with stanford data? ----
barcodes <- get_biodata("ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/CSF_plasma_matched_barcodes.csv")
CSF_patients <- get_biodata("ADRC/ADRC_Data_Sharing_2023/CSF_metadata_2023-03-04.csv")
CSF_prots <- get_biodata("ADRC/CSF_MedNorm/input_DE_HumanSoma_SamplesOnly_LOD/CSFProts.log10.noLODFilter.csv")
plasma_patients <- get_biodata("ADRC/ADRC_Data_Sharing_2023/Plasma_metadata_samplesOnly_2023-03-04.csv")
plasma_prots <- get_biodata("ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/plasmaProts.log10.csv")

CSF_patients$Age <- gsub('^"&"$', '', CSF_patients$Age)
CSF_patients$Age <- as.numeric(CSF_patients$Age)

common_prots <- intersect(names(CSF_prots), names(plasma_prots))
CSF_prots <- data.frame(Barcode = CSF_patients$Barcode, CSF_prots[, common_prots]) 
plasma_prots <- data.frame(Barcode = plasma_prots$Barcode, plasma_prots[, common_prots])

CSF_prots_fil <- CSF_prots[CSF_prots$Barcode %in% barcodes$CSF_Barcode, ]
CSF_patients_fil <- CSF_patients[CSF_patients$Barcode %in% barcodes$CSF_Barcode, ]
plasma_prots_fil <- plasma_prots %>%
  filter(Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, barcodes$Plasma_Barcode)))

age_fil <- plasma_patients %>% 
  dplyr::select(Age, Barcode) %>% 
  rename(Plasma_Barcode = Barcode) %>% 
  filter(Plasma_Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes$Plasma_Barcode)) %>% 
  mutate(CSF_Barcode = CSF_patients_fil$Barcode) %>% 
  dplyr::select(Age, CSF_Barcode, Plasma_Barcode)

patientdata_fil <- plasma_patients %>% 
  rename(Plasma_Barcode = Barcode, Plasma_Zscore = ConnectivityZscore) %>% 
  filter(Plasma_Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes$Plasma_Barcode)) 

gendermatch_index <- CSF_patients_fil$Gender == patientdata_fil$Gender
CSF_patients_fil <- CSF_patients_fil[gendermatch_index, ]
patientdata_fil <- patientdata_fil[gendermatch_index, ]
age_fil <- age_fil[gendermatch_index, ]
CSF_prots_fil <- CSF_prots_fil[gendermatch_index, ]
plasma_prots_fil <- plasma_prots_fil[gendermatch_index, ]

patientdata_fil <- patientdata_fil %>% 
  mutate(CSF_Barcode = CSF_patients_fil$Barcode, Plasma_drawage = Age, CSF_drawage = CSF_patients_fil$Age, 
         Plasma_storagedays = Storage_days, CSF_storagedays = CSF_patients_fil$Storage_days,
         CSF_days = CSF_patients_fil$Storage_days, Plasma_drawdate = strptime(Date.of.draw, format = "%m/%d/%y"), 
         CSF_drawdate = strptime(CSF_patients_fil$Date_of_CSF_sample, format = "%m/%d/%y"), 
         Plasma_drawstatus = Diagnosis_group, CSF_drawstatus = CSF_patients_fil$Diagnosis_group) %>% 
  group_by(Plasma_Barcode) %>%
  mutate(avg_drawage = (Plasma_drawage + as.numeric(CSF_drawage)) / 2, 
         drawdate_diff = difftime(Plasma_drawdate, CSF_drawdate, units = "days"), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>%
  dplyr::select(Plasma_Barcode, CSF_Barcode, Gender, avg_drawage, Plasma_drawdate, Plasma_drawage, Plasma_drawstatus, 
                Plasma_storagedays, CSF_drawdate, CSF_drawage, CSF_drawstatus, CSF_storagedays, drawdate_diff, avg_drawdate) %>%
  as.data.frame()

# actually drawdate diff is all below 90!!! we can proceed without setting a time window. 
patientdata_fil <- patientdata_fil %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "ADMCI", "AD", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "MCI-AD", "AD", CSF_drawstatus)) %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "HC", "CO", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "HC", "CO", CSF_drawstatus)) %>% 
  mutate(final_status = 
           if_else(Plasma_drawdate < CSF_drawdate, Plasma_drawstatus, CSF_drawstatus))

patientdata_AD_index <- patientdata_fil$final_status == "AD" 
patientdata_CO_index <- patientdata_fil$final_status == "CO"
patientdata_AD_index[is.na(patientdata_AD_index)] <- FALSE
patientdata_CO_index[is.na(patientdata_CO_index)] <- FALSE
patientdata_index <- patientdata_AD_index | patientdata_CO_index

# clean WashU data 
plasma_meta_temp <- get_biodata("WashU_ADRC/Plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
plasma_prots_temp <- get_biodata("WashU_ADRC/Plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
plasma_patient_temp <- read.xlsx(paste(dir, "WashU_ADRC/Plasma/Plasma_ BasicDemo.xlsx", sep = "/"))
CSF_meta_temp <- get_biodata("WashU_ADRC/CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
CSF_prots_temp <- get_biodata("WashU_ADRC/CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
CSF_patient_temp <- read.xlsx(paste(dir, "WashU_ADRC/CSF/CSF_BasicDemo.xlsx", sep = "/"))

patientmeta_temp <- get_biodata("WashU_ADRC/commonfile/WashU_ADNI_demographic.csv")

plasma_prots_temp[4:ncol(plasma_prots_temp)] <- log10(plasma_prots_temp[4:ncol(plasma_prots_temp)])
CSF_prots_temp[4:ncol(CSF_prots_temp)] <- log10(CSF_prots_temp[4:ncol(CSF_prots_temp)])

protMeta <- read_csv("/labs/twc/jarod/Data/ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/ProteinMetadata_with_LOD.csv", col_types = )
protMeta$SeqId = str_replace(protMeta$SeqId, "-", ".")
protMeta <- filter(protMeta, (Organism == "Human" | Organism == "HIV-1" | Organism == "HIV-2") & Type == "Protein")
all_prots_ID <- names(plasma_prots_temp)[4:ncol(plasma_prots_temp)] 
all_prots <- str_remove(all_prots_ID, "^X")
my_order <- match(protMeta$SeqId, all_prots)
all_prots <- protMeta$Key_2
nprots <- length(all_prots)

plasma_prots_temp <- plasma_prots_temp[c(1:3, 3 + my_order)]
colnames(plasma_prots_temp) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)

CSF_prots_temp <- CSF_prots_temp[c(1:3, 3 + my_order)]
colnames(CSF_prots_temp) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)


all_prots <- intersect(names(plasma_prots_temp)[4:ncol(plasma_prots_temp)], names(plasma_prots)[4:ncol(plasma_prots)])
nprots <- length(all_prots)

plasma_patient_temp <- plasma_patient_temp %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed AD", 
                   "AD", final_cc_status.updated)) %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed Control", 
                   "CO", final_cc_status.updated)) %>% 
  mutate(DateOfBirth = if_else(grepl("/", DateOfBirth), strptime(DateOfBirth, format = "%m/%d/%Y"), 
                               convertToDateTime(as.numeric(DateOfBirth)))) %>% 
  mutate(drawdate = convertToDateTime(as.numeric(drawdate)))

CSF_patient_temp <- CSF_patient_temp %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed AD", 
                   "AD", Last.status)) %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed Control", 
                   "CO", Last.status)) %>% 
  mutate(DrawDate = convertToDateTime(DrawDate))

plasma_meta_temp <- arrange(plasma_meta_temp, PA_DB_UID)
plasma_prots_temp <- arrange(plasma_prots_temp, PA_DB_UID)
plasma_patient_temp <- arrange(plasma_patient_temp, PA_DB_UID)
CSF_meta_temp <- arrange(CSF_meta_temp, PA_DB_UID)
CSF_prots_temp <- arrange(CSF_prots_temp, PA_DB_UID)
CSF_patient_temp <- arrange(CSF_patient_temp, PA_DB_UID)

common_ID <- intersect(CSF_patient_temp$PA_DB_UID, plasma_patient_temp$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .) %>% arrange(PA_DB_UID) 
common_ID <- common_ID$PA_DB_UID

plasma_meta_temp <- plasma_meta_temp[plasma_meta_temp$PA_DB_UID %in% common_ID, ]
plasma_prots_temp <- plasma_prots_temp[plasma_prots_temp$PA_DB_UID %in% common_ID, ]
plasma_patient_temp <- plasma_patient_temp[plasma_patient_temp$PA_DB_UID %in% common_ID, ]
CSF_meta_temp <- CSF_meta_temp[CSF_meta_temp$PA_DB_UID %in% common_ID, ]
CSF_prots_temp <- CSF_prots_temp[CSF_prots_temp$PA_DB_UID %in% common_ID, ]
CSF_patient_temp <- CSF_patient_temp[CSF_patient_temp$PA_DB_UID %in% common_ID, ]

patientdata_temp_fil <- merge(CSF_patient_temp, plasma_patient_temp, by = "PA_DB_UID") %>% 
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

date_window <- 120
date_index <- patientdata_temp_fil$drawdate_diff <= date_window
CSF_prots_temp_fil <- CSF_prots_temp[date_index, ]
plasma_prots_temp_fil <- plasma_prots_temp[date_index, ]

patientdata_temp <- patientdata_temp_fil[date_index, ] %>% 
  mutate(final_status = ifelse(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus), 
         avg_drawage = abs(Plasma_drawage + CSF_drawage) / 2) %>% 
  select(PA_DB_UID, Sex, avg_drawage, drawdate_diff, Plasma_drawdate, CSF_drawdate, final_status, storage_days) 


# get rid of non AD or COs
knight_index <- patientdata_temp$final_status == "AD" | patientdata_temp$final_status == "CO" 
stanford_index <- patientdata_fil$final_status == "AD" | patientdata_fil$final_status == "CO"

patientdata_temp <- patientdata_temp[knight_index, ]
patientdata_fil <- patientdata_fil[stanford_index, ]
plasma_prots_temp_fil <- plasma_prots_temp_fil[knight_index, ]
plasma_prots_fil <- plasma_prots_fil[stanford_index, ]
CSF_prots_temp_fil <- CSF_prots_temp_fil[knight_index, ]
CSF_prots_fil <- CSF_prots_fil[stanford_index, ]

all_patientdata <- bind_rows(patientdata_fil %>% mutate(Sex = Gender, drawdate_diff = abs(drawdate_diff)) %>% 
                               dplyr::select(Sex, avg_drawage, final_status, drawdate_diff), 
                             patientdata_temp %>% 
                               dplyr::select(Sex, avg_drawage, final_status, drawdate_diff))
all_CSF <- bind_rows(CSF_prots_fil[-c(1:3)], CSF_prots_temp_fil[all_prots])
all_plasma <- bind_rows(plasma_prots_fil[-c(1:3)], plasma_prots_temp_fil[all_prots])

all_CSF <- all_CSF[!is.na(all_patientdata$Sex), ]
all_plasma <- all_plasma[!is.na(all_patientdata$Sex), ]
all_patientdata <- all_patientdata[!is.na(all_patientdata$Sex), ]

AD_index <- all_patientdata$final_status == "AD"
CO_index <- all_patientdata$final_status == "CO"

all_patientdata_AD <- all_patientdata[AD_index, ]
all_patientdata_CO <- all_patientdata[CO_index, ]

all_ratio <- all_CSF / all_plasma
all_ratio_AD <- all_ratio[all_patientdata$final_status == "AD", ]
all_ratio_CO <- all_ratio[all_patientdata$final_status == "CO", ]

all_ratio_z <- scale(all_ratio)
all_ratio_AD_z <- scale(all_ratio_AD)
all_ratio_CO_z <- scale(all_ratio_CO)

albumin_all_ratio <- all_ratio[, names(all_ratio) == "ALB.18380.78.3"]
albumin_all_ratio_z <- scale(albumin_all_ratio)


# Create histogram of population 
drawdate_windows <- 20 * (0:6)
generate_sampleHist(all_patientdata, drawdate_windows) 

# predict data and trends 
age_min <- max(min(all_patientdata_AD$avg_drawage), min(all_patientdata_CO$avg_drawage))
age_max <- min(max(all_patientdata_AD$avg_drawage), max(all_patientdata_CO$avg_drawage))
age_seq <- seq(age_min, age_max, by = 0.25)

predratio_AD_z <- loess_predict(all_patientdata_AD$avg_drawage, all_ratio_AD_z, age_seq, "Age")
predratio_CO_z <- loess_predict(all_patientdata_CO$avg_drawage, all_ratio_CO_z, age_seq, "Age")

# clustering for CO and DA
ratio_CO_z_dist <- dist(t(predratio_CO_z[-1]), method = "euclidean")
ratio_CO_z_clust <- hclust(ratio_CO_z_dist, method = "complete")
ratio_AD_z_dist <- dist(t(predratio_AD_z[-1]), method = "euclidean")
ratio_AD_z_clust <- hclust(ratio_AD_z_dist, method = "complete")

determine_numclust(ratio_CO_z_clust, ratio_CO_z_dist)
determine_numclust(ratio_AD_z_clust, ratio_AD_z_dist)
# we conclude both 9 clusters 

nclust <- 10

ratio_AD_z_clustnum <- proteinByClust(ratio_AD_z_clust, ratio_AD_z_dist, nclust)
ratio_CO_z_clustnum <- proteinByClust(ratio_CO_z_clust, ratio_CO_z_dist, nclust) # tbh 8 clusters are better

ratio_AD_z_long <- proteinClustData(predratio_AD_z, ratio_AD_z_clustnum)
ratio_CO_z_long <- proteinClustData(predratio_CO_z, ratio_CO_z_clustnum)

generate_clusterplot(ratio_AD_z_long)
generate_clusterplot(ratio_CO_z_long)
ratio_CO_z_clustnum[ratio_CO_z_clustnum$cluster == 10, ]


# LOOK AT RATIO BW PLASMA AND CSF ALBUMIN (albumin index). look at which ones 
# are correlated with albumin. (ALB.18380.78.3) Find which cluster it is in, etc 
albumin_cluster_AD <- ratio_AD_z_clustnum[ratio_AD_z_clustnum$protein == "ALB.18380.78.3", ]$cluster
albumin_prots_AD <-ratio_AD_z_clustnum[ratio_AD_z_clustnum$cluster == albumin_cluster_AD, ]$protein
albumin_cluster_CO <- ratio_CO_z_clustnum[ratio_CO_z_clustnum$protein == "ALB.18380.78.3", ]$cluster
albumin_prots_CO <-ratio_CO_z_clustnum[ratio_CO_z_clustnum$cluster == albumin_cluster_CO, ]$protein
albumin_intersection <- intersect(albumin_prots_AD, albumin_prots_CO)

# stratified analysis: remake the volcano plots 
ratio_AD_models <- generate_lmodels(all_ratio_AD, all_patientdata_AD)
ratio_CO_models <- generate_lmodels(all_ratio_CO, all_patientdata_CO)

ratio_AD_summary <- generate_lmsummary(ratio_AD_models)
ratio_CO_summary <- generate_lmsummary(ratio_CO_models)

xVars <- c("Male", "Age")

ratio_AD_volcanodata <- generate_volcanodata(ratio_AD_summary)
ratio_CO_volcanodata <- generate_volcanodata(ratio_CO_summary)


ratio_AD_volcanodata_nonpadj <- generate_volcanodata_nonpadj(ratio_AD_summary)

generate_volcanoplot(ratio_AD_volcanodata)
generate_volcanoplot(ratio_CO_volcanodata)
generate_volcanoplot(ratio_AD_volcanodata_nonpadj)  


# Get enrichR pathways 
albumin_enriched <- enrichr(str_extract(albumin_intersection, "\\w+(?=\\.)"), "ChEA_2022")
albumin_enriched <- albumin_enriched$ChEA_2022
albumin_AD_enriched <- enrichr(str_extract(albumin_prots_AD, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022
albumin_CO_enriched <- enrichr(str_extract(albumin_prots_CO, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022

ratio_CO_sex_up <- getTopProteins(ratio_CO_volcanodata[[1]], "up")
ratio_CO_sex_down <- getTopProteins(ratio_CO_volcanodata[[1]], "down")
ratio_CO_age_up <- getTopProteins(ratio_CO_volcanodata[[2]], "up")
ratio_CO_age_down <- getTopProteins(ratio_CO_volcanodata[[2]], "down")

ratio_AD_sex_up <- getTopProteins(ratio_AD_volcanodata[[1]], "up")
ratio_AD_sex_down <- getTopProteins(ratio_AD_volcanodata[[1]], "down")
ratio_AD_age_up <- getTopProteins(ratio_AD_volcanodata[[2]], "up")
ratio_AD_age_down <- getTopProteins(ratio_AD_volcanodata[[2]], "down")

ratio_CO_sex_up_enriched <- enrichVolcanodata(ratio_CO_sex_up, "ChEA_2022")$ChEA_2022
ratio_CO_sex_down_enriched <- enrichVolcanodata(ratio_CO_sex_down, "ChEA_2022")$ChEA_2022
ratio_CO_age_up_enriched <- enrichVolcanodata(ratio_CO_age_up, "ChEA_2022")$ChEA_2022
ratio_CO_age_down_enriched <- enrichVolcanodata(ratio_CO_age_down, "ChEA_2022")$ChEA_2022

ratio_AD_sex_up_enriched <- enrichVolcanodata(ratio_AD_sex_up, "ChEA_2022")$ChEA_2022
ratio_AD_sex_down_enriched <- enrichVolcanodata(ratio_AD_sex_down, "ChEA_2022")$ChEA_2022
ratio_AD_age_up_enriched <- enrichVolcanodata(ratio_AD_age_up, "ChEA_2022")$ChEA_2022
ratio_AD_age_down_enriched <- enrichVolcanodata(ratio_AD_age_down, "ChEA_2022")$ChEA_2022


# cluster preservation analysis 
## centroid matrix: row = each gene, col = each cluster


# knight before p-adjusted VS stanford volcano plots --> compare enrichR pathways 
knight_ratio <- CSF_prots_temp_fil[-c(1:3)] / plasma_prots_temp_fil[-c(1:3)]
knight_all_prots <- colnames(knight_ratio)
knight_patients <- patientdata_temp

knight_AD_ratio <- knight_ratio[knight_patients$final_status == "AD", ]
knight_CO_ratio <- knight_ratio[knight_patients$final_status == "CO", ]
knight_AD_patients <- knight_patients[knight_patients$final_status == "AD", ]
knight_CO_patients <- knight_patients[knight_patients$final_status == "CO", ]

knight_AD_lmodels <- generate_lmodels(knight_AD_ratio, knight_AD_patients, c("Sex", "avg_drawage"))
knight_CO_lmodels <- generate_lmodels(knight_CO_ratio, knight_CO_patients, c("Sex", "avg_drawage"))

knight_AD_lmsummary <- generate_lmsummary(knight_AD_lmodels, knight_all_prots, c("Male", "Age"))
knight_CO_lmsummary <- generate_lmsummary(knight_CO_lmodels, knight_all_prots, c("Male", "Age"))

knight_AD_volcanodata <- generate_volcanodata(knight_AD_lmsummary, c("Male", "Age"))
knight_CO_volcanodata <- generate_volcanodata(knight_CO_lmsummary, c("Male", "Age"))

generate_volcanoplot(knight_AD_volcanodata)
generate_volcanoplot(knight_CO_volcanodata)

#interaction model bw age and alzheimers disease status (extract CD? score)
knight_patients_CDR <- plasma_patient_temp[plasma_patient_temp$PA_DB_UID %in% knight_patients$PA_DB_UID, ] %>%
  select(CDR_at_blood_draw)
knight_patients_CDR <- knight_patients %>% mutate(CDR = knight_patients_CDR$CDR_at_blood_draw) %>%
  mutate(Age_CDR = avg_drawage * CDR)

knight_interaction_lmodels <- generate_lmodels(knight_ratio, knight_patients_CDR, c("Sex", "Age_CDR"))
knight_interaction_lmsummary <- generate_lmsummary(knight_interaction_lmodels, knight_all_prots, c("Male", "Age*CDR"))
knight_interaction_volcanodata <- generate_volcanodata(knight_interaction_lmsummary, c("Male", "Age*CDR"))
generate_volcanoplot(knight_interaction_volcanodata, c("Male", "Age*CDR"))

## compare with combined dataset
all_interaction_lmodels <- list()
for(i in 1:ncol(all_ratio)) {
  all_interaction_lmodels[[i]] <- summary(lm(all_ratio[, i] ~ all_patientdata$Sex 
                                             + all_patientdata$avg_drawage 
                                             + relevel(factor(all_patientdata$final_status), ref = "CO")))
}
all_prots <- names(all_ratio)
all_interaction_lmsummary <- generate_lmsummary(all_interaction_lmodels, all_prots, c("Male", "Age", "AD"))
all_interaction_volcanodata <- generate_volcanodata(all_interaction_lmsummary, c("Male", "Age", "AD"))
generate_volcanoplot(all_interaction_volcanodata, c("Male", "Age", "AD"))

## enrichment analysis 
all_ADCO_sex_up <- getTopProteins(all_interaction_volcanodata[[1]], "up")
all_ADCO_sex_down <- getTopProteins(all_interaction_volcanodata[[1]], "down")
all_ADCO_age_up <- getTopProteins(all_interaction_volcanodata[[2]], "up")
all_ADCO_age_down <- getTopProteins(all_interaction_volcanodata[[2]], "down")
all_ADCO_AD_up <- getTopProteins(all_interaction_volcanodata[[3]], "up")
all_ADCO_AD_down <- getTopProteins(all_interaction_volcanodata[[3]], "down")

all_ADCO_sex_up_enriched <- enrichVolcanodata(all_ADCO_sex_up, "ChEA_2022")$ChEA_2022
all_ADCO_sex_down_enriched <- enrichVolcanodata(all_ADCO_sex_down, "ChEA_2022")$ChEA_2022
all_ADCO_age_up_enriched <- enrichVolcanodata(all_ADCO_age_up, "ChEA_2022")$ChEA_2022
all_ADCO_age_down_enriched <- enrichVolcanodata(all_ADCO_age_down, "ChEA_2022")$ChEA_2022
all_ADCO_AD_up_enriched <- enrichVolcanodata(all_ADCO_AD_up, "ChEA_2022")$ChEA_2022
all_ADCO_AD_down_enriched <- enrichVolcanodata(all_ADCO_AD_down, "ChEA_2022")$ChEA_2022

knight_ADCO_sex_up <- getTopProteins(knight_interaction_volcanodata[[1]], "up")
knight_ADCO_sex_down <- getTopProteins(knight_interaction_volcanodata[[1]], "down")
knight_ADCO_int_up <- getTopProteins(knight_interaction_volcanodata[[2]], "up")
knight_ADCO_int_down <- getTopProteins(knight_interaction_volcanodata[[2]], "down")

knight_ADCO_sex_up_enriched <- enrichVolcanodata(knight_ADCO_sex_up, "ChEA_2022")$ChEA_2022
knight_ADCO_sex_down_enriched <- enrichVolcanodata(knight_ADCO_sex_down, "ChEA_2022")$ChEA_2022
knight_ADCO_int_up_enriched <- enrichVolcanodata(knight_ADCO_int_up, "ChEA_2022")$ChEA_2022
knight_ADCO_int_down_enriched <- enrichVolcanodata(knight_ADCO_int_down, "ChEA_2022")$ChEA_2022

# Clustering and enrichment analysis 
ratio_AD_clustable <- table(ratio_AD_z_clustnum$cluster)
ratio_CO_clustable <- table(ratio_CO_z_clustnum$cluster)

ratio_AD_clust_enriched <- enrichClusts(ratio_AD_z_clustnum)
ratio_CO_clust_enriched <- enrichClusts(ratio_CO_z_clustnum)


# differential sliding window analysis
sliding_window <- 10

