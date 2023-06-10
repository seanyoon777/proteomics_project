rm(list = ls())

source("Scripts/Helpers/init.R")
source("Scripts/Helpers/differentialAnalysis.R")
source("Scripts/Helpers/hierClust.R")
source("Scripts/Helpers/demoHist.R")
source("Scripts/Helpers/enrichAnalysis.R")

dir <- "/labs/twc/jarod/Data"
get_biodata <- function(path) {
  read.csv(paste(dir, path, sep = "/"))
}


# same thing with stanford data? ----
stanford_barcodes <- get_biodata("ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/CSF_plasma_matched_barcodes.csv")
stanford_CSF_patients <- get_biodata("ADRC/ADRC_Data_Sharing_2023/CSF_metadata_2023-03-04.csv")
stanford_CSF_prots <- get_biodata("ADRC/CSF_MedNorm/input_DE_HumanSoma_SamplesOnly_LOD/CSFProts.log10.noLODFilter.csv")
stanford_plasma_patients <- get_biodata("ADRC/ADRC_Data_Sharing_2023/Plasma_metadata_samplesOnly_2023-03-04.csv")
stanford_plasma_prots <- get_biodata("ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/plasmaProts.log10.csv")

stanford_CSF_patients$Age <- gsub('^"&"$', '', stanford_CSF_patients$Age)
stanford_CSF_patients$Age <- as.numeric(stanford_CSF_patients$Age)

stanford_plasma_prots <- stanford_plasma_prots[stanford_plasma_prots$Barcode %in% stanford_plasma_patients$Barcode, ]
stanford_CSF_prots <- data.frame(Barcode = stanford_CSF_patients$Barcode, stanford_CSF_prots)

stanford_common_prots <- intersect(names(stanford_CSF_prots), names(stanford_plasma_prots))

stanford_plasma_prots <- stanford_plasma_prots[stanford_common_prots]
stanford_CSF_prots <- stanford_CSF_prots[stanford_common_prots]

stanford_plasma_prots <- stanford_plasma_prots[stanford_plasma_prots$Barcode %in% stanford_barcodes$Plasma_Barcode, ]
stanford_plasma_patients <- stanford_plasma_patients[stanford_plasma_patients$Barcode %in% stanford_barcodes$Plasma_Barcode, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_CSF_prots$Barcode %in% stanford_barcodes$CSF_Barcode, ]
stanford_CSF_patients <- stanford_CSF_patients[stanford_CSF_patients$Barcode %in% stanford_barcodes$CSF_Barcode, ]

stanford_plasma_prots <- stanford_plasma_prots %>%
  filter(Barcode %in% stanford_barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, stanford_barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, stanford_barcodes$Plasma_Barcode)))

stanford_plasma_patients <- stanford_plasma_patients %>%
  filter(Barcode %in% stanford_barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, stanford_barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, stanford_barcodes$Plasma_Barcode)))

stanford_gendermatch_index <- stanford_CSF_patients$Gender == stanford_plasma_patients$Gender
stanford_CSF_patients <- stanford_CSF_patients[stanford_gendermatch_index, ]
stanford_plasma_patients <- stanford_plasma_patients[stanford_gendermatch_index, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_gendermatch_index, ]
stanford_plasma_prots <- stanford_plasma_prots[stanford_gendermatch_index, ]

stanford_patients <- stanford_plasma_patients %>% 
  mutate(Plasma_Barcode = Barcode, CSF_Barcode = stanford_CSF_patients$Barcode, 
         Plasma_drawage = Age, CSF_drawage = stanford_CSF_patients$Age, 
         Plasma_storagedays = Storage_days, 
         CSF_storagedays = stanford_CSF_patients$Storage_days,
         CSF_days = stanford_CSF_patients$Storage_days, 
         Plasma_drawdate = strptime(Date.of.draw, format = "%m/%d/%y"), 
         CSF_drawdate = strptime(stanford_CSF_patients$Date_of_CSF_sample, format = "%m/%d/%y"), 
         Plasma_drawstatus = Diagnosis_group, 
         CSF_drawstatus = stanford_CSF_patients$Diagnosis_group) %>% 
  group_by(Plasma_Barcode) %>%
  mutate(avg_drawage = (Plasma_drawage + as.numeric(CSF_drawage)) / 2, 
         drawdate_diff = difftime(Plasma_drawdate, CSF_drawdate, units = "days"), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>%
  mutate(Plasma_drawdate2 = as.numeric(difftime(today(), Plasma_drawdate)), batch_effect = 1) %>%
  dplyr::select(Plasma_Barcode, CSF_Barcode, Gender, avg_drawage, Plasma_drawdate, 
                Plasma_drawage, Plasma_drawstatus, Plasma_storagedays, CSF_drawdate, 
                CSF_drawage, CSF_drawstatus, CSF_storagedays, drawdate_diff,
                Plasma_drawdate2, avg_drawdate, batch_effect) %>%
  as.data.frame()

# actually drawdate diff is all below 90!!! we can proceed without setting a time window. 
stanford_patients <- stanford_patients %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "ADMCI", "AD", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "MCI-AD", "AD", CSF_drawstatus)) %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "HC", "CO", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "HC", "CO", CSF_drawstatus)) %>% 
  mutate(final_status = 
           if_else(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus))

stanford_AD_index <- stanford_patients$final_status == "AD" 
stanford_CO_index <- stanford_patients$final_status == "CO"
stanford_AD_index[is.na(stanford_AD_index)] <- FALSE
stanford_CO_index[is.na(stanford_CO_index)] <- FALSE
stanford_index <- stanford_AD_index | stanford_CO_index

# clean WashU data 
knight_plasma_meta <- get_biodata("WashU_ADRC/Plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
knight_plasma_prots <- get_biodata("WashU_ADRC/Plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
knight_plasma_patients <- read.xlsx(paste(dir, "WashU_ADRC/Plasma/Plasma_ BasicDemo.xlsx", sep = "/"))
knight_CSF_meta <- get_biodata("WashU_ADRC/CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
knight_CSF_prots <- get_biodata("WashU_ADRC/CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
knight_CSF_patients <- read.xlsx(paste(dir, "WashU_ADRC/CSF/CSF_BasicDemo.xlsx", sep = "/"))

knight_patients_meta <- get_biodata("WashU_ADRC/commonfile/WashU_ADNI_demographic.csv")

knight_plasma_prots[4:ncol(knight_plasma_prots)] <- log10(knight_plasma_prots[4:ncol(knight_plasma_prots)])
knight_CSF_prots[4:ncol(knight_CSF_prots)] <- log10(knight_CSF_prots[4:ncol(knight_CSF_prots)])

protMeta <- read_csv("/labs/twc/jarod/Data/ADRC/Plasma_ANML/input_DE_HumanSomas_Samples_Only_LOD/ProteinMetadata_with_LOD.csv", col_types = )
protMeta$SeqId = str_replace(protMeta$SeqId, "-", ".")
protMeta <- filter(protMeta, (Organism == "Human" | Organism == "HIV-1" | Organism == "HIV-2") & Type == "Protein")
all_prots_ID <- names(knight_plasma_prots)[4:ncol(knight_plasma_prots)] 
all_prots <- str_remove(all_prots_ID, "^X")
my_order <- match(protMeta$SeqId, all_prots)
all_prots <- protMeta$Key_2
nprots <- length(all_prots)

knight_plasma_prots <- knight_plasma_prots[c(1:3, 3 + my_order)]
colnames(knight_plasma_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)

knight_CSF_prots <- knight_CSF_prots[c(1:3, 3 + my_order)]
colnames(knight_CSF_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)


all_prots <- intersect(names(knight_plasma_prots)[4:ncol(knight_plasma_prots)], 
                       names(knight_plasma_prots)[4:ncol(knight_plasma_prots)])
nprots <- length(all_prots)

knight_plasma_patients <- knight_plasma_patients %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed AD", 
                   "AD", final_cc_status.updated)) %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed Control", 
                   "CO", final_cc_status.updated)) %>% 
  mutate(DateOfBirth = if_else(grepl("/", DateOfBirth), strptime(DateOfBirth, format = "%m/%d/%Y"), 
                               convertToDateTime(as.numeric(DateOfBirth)))) %>% 
  mutate(drawdate = convertToDateTime(as.numeric(drawdate)))

knight_CSF_patients <- knight_CSF_patients %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed AD", 
                   "AD", Last.status)) %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed Control", 
                   "CO", Last.status)) %>% 
  mutate(DrawDate = convertToDateTime(DrawDate))

knight_plasma_meta <- arrange(knight_plasma_meta, PA_DB_UID)
knight_plasma_prots <- arrange(knight_plasma_prots, PA_DB_UID)
knight_plasma_patients <- arrange(knight_plasma_patients, PA_DB_UID)
knight_CSF_meta <- arrange(knight_CSF_meta, PA_DB_UID)
knight_CSF_prots <- arrange(knight_CSF_prots, PA_DB_UID)
knight_CSF_patients <- arrange(knight_CSF_patients, PA_DB_UID)

knight_common_ID <- intersect(knight_CSF_patients$PA_DB_UID, knight_plasma_patients$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .) %>% arrange(PA_DB_UID) 
knight_common_ID <- knight_common_ID$PA_DB_UID

knight_plasma_meta <- knight_plasma_meta[knight_plasma_meta$PA_DB_UID %in% knight_common_ID, ]
knight_plasma_prots <- knight_plasma_prots[knight_plasma_prots$PA_DB_UID %in% knight_common_ID, ]
knight_plasma_patients <- knight_plasma_patients[knight_plasma_patients$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_meta <- knight_CSF_meta[knight_CSF_meta$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_prots <- knight_CSF_prots[knight_CSF_prots$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_patients <- knight_CSF_patients[knight_CSF_patients$PA_DB_UID %in% knight_common_ID, ]

knight_patients <- merge(knight_CSF_patients, knight_plasma_patients, by = "PA_DB_UID") %>% 
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status) %>%
  group_by(PA_DB_UID) %>% 
  mutate(drawdate_diff = abs(difftime(Plasma_drawdate, CSF_drawdate, units = "days")), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>% 
  mutate(storage_days = as.numeric(difftime(today(), Plasma_drawdate, units = "days")), batch_effect = 0) %>% 
  as.data.frame() 

date_window <- 120
knight_date_index <- knight_patients$drawdate_diff <= date_window
knight_CSF_prots <- knight_CSF_prots[knight_date_index, ]
knight_plasma_prots <- knight_plasma_prots[knight_date_index, ]

knight_patients <- knight_patients[knight_date_index, ] %>% 
  mutate(final_status = ifelse(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus), 
         avg_drawage = abs(Plasma_drawage + CSF_drawage) / 2, storage_days2 = storage_days) %>% 
  select(PA_DB_UID, Sex, avg_drawage, drawdate_diff, Plasma_drawdate, CSF_drawdate, 
         final_status, storage_days, storage_days2, batch_effect) 


# get rid of non AD or COs
knight_index <- knight_patients$final_status == "AD" | knight_patients$final_status == "CO" 
stanford_index <- stanford_patients$final_status == "AD" | stanford_patients$final_status == "CO"

knight_patients <- knight_patients[knight_index, ]
stanford_patients <- stanford_patients[stanford_index, ]
knight_plasma_prots <- knight_plasma_prots[knight_index, ]
stanford_plasma_prots <- stanford_plasma_prots[stanford_index, ]
knight_CSF_prots <- knight_CSF_prots[knight_index, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_index, ]

all_patientdata <- bind_rows(stanford_patients %>% 
                               mutate(Sex = Gender, drawdate_diff = abs(drawdate_diff)) %>% 
                               mutate(storage_days = Plasma_drawdate2) %>% 
                               dplyr::select(Sex, avg_drawage, final_status, drawdate_diff, 
                                             storage_days, batch_effect), 
                             knight_patients %>% 
                               dplyr::select(Sex, avg_drawage, final_status, drawdate_diff, 
                                             storage_days, batch_effect))

all_prots <- intersect(names(stanford_CSF_prots), names(knight_CSF_prots))
all_CSF <- bind_rows(stanford_CSF_prots[all_prots], knight_CSF_prots[all_prots])
all_plasma <- bind_rows(stanford_plasma_prots[all_prots], knight_plasma_prots[all_prots])

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
xVars <- c("Sex", "avg_drawage", "drawdate_diff", "storage_days", "batch_effect")
ratio_AD_models <- generate_lmodels(all_ratio_AD, all_patientdata_AD, xVars)
ratio_CO_models <- generate_lmodels(all_ratio_CO, all_patientdata_CO, xVars)

xLabs <- c("Male", "Age", "Drawdate difference", "Storage days", "Study bias")
ratio_AD_summary <- generate_lmsummary(ratio_AD_models, colnames(all_ratio), xLabs)
ratio_CO_summary <- generate_lmsummary(ratio_CO_models, colnames(all_ratio), xLabs)

ratio_AD_volcanodata <- generate_volcanodata(ratio_AD_summary, xLabs)
ratio_CO_volcanodata <- generate_volcanodata(ratio_CO_summary, xLabs)


ratio_AD_volcanodata_nonpadj <- generate_volcanodata_nonpadj(ratio_AD_summary)

xLabs <- c("Male", "Age")
generate_volcanoplot(ratio_AD_volcanodata, xLabs, 3)
generate_volcanoplot(ratio_CO_volcanodata, xLabs, 3)
generate_volcanoplot(ratio_AD_volcanodata_nonpadj, xLabs, 3)

xVars <- c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect")
ratio_models <- list()
for(i in 1:ncol(all_ratio)) {
  ratio_models[[i]] <- summary(lm(all_ratio[, i] ~ all_patientdata$Sex + all_patientdata$avg_drawage + 
                               relevel(factor(all_patientdata$final_status), ref = "CO") + 
                               all_patientdata$drawdate_diff + 
                               all_patientdata$storage_days + 
                               all_patientdata$batch_effect))
}
xLabs <- c("Male", "Age", "AD", "Drawdate difference", "storage days", "Study bias")
ratio_summary <- generate_lmsummary(ratio_models, colnames(all_ratio), xLabs)
ratio_volcanodata <- generate_volcanodata(ratio_summary, xLabs)
generate_volcanoplot(ratio_volcanodata, c("Male", "Age", "AD"), 3)

xVars <- c("sex", "age", "AD", "drawdateDiff", "storageDays", "studyBias")
for(i in 1:6) {
  upProteins <- getTopProteins(ratio_volcanodata[[i]], "up")
  downProteins <- getTopProteins(ratio_volcanodata[[i]], "down")
  upEnriched <- enrichVolcanodata(upProteins, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_ADCOCombined/combinedData_", xVars[i], "_upRegulated.csv"))
  downEnriched <- enrichVolcanodata(downProteins, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_ADCOCombined/combinedData_", xVars[i], "_downRegulated.csv"))
}

# Get enrichR pathways 
albumin_enriched <- enrichr(str_extract(albumin_intersection, "\\w+(?=\\.)"), "ChEA_2022")
albumin_enriched <- albumin_enriched$ChEA_2022
albumin_AD_enriched <- enrichr(str_extract(albumin_prots_AD, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022
albumin_CO_enriched <- enrichr(str_extract(albumin_prots_CO, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022
write.csv(albumin_AD_enriched, paste0(getwd(), "/Data/Cluster/albumin_AD.csv"))
write.csv(albumin_CO_enriched, paste0(getwd(), "/Data/Cluster/albumin_CO.csv"))

xVars <- c("sex", "age", "drawdateDiff", "storageDays", "studyBias")
for(i in 1:5) {
  upProteinsCO <- getTopProteins(ratio_CO_volcanodata[[i]], "up")
  downProteinsCO <- getTopProteins(ratio_CO_volcanodata[[i]], "down")
  upEnriched <- enrichVolcanodata(upProteinsCO, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_CO/combinedData_CO_", xVars[i], "_upRegulated.csv"))
  downEnriched <- enrichVolcanodata(downProteinsCO, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_CO/combinedData_CO_", xVars[i], "_downRegulated.csv"))
  
  upProteinsAD <- getTopProteins(ratio_AD_volcanodata[[i]], "up")
  downProteinsAD <- getTopProteins(ratio_AD_volcanodata[[i]], "down")
  upEnriched <- enrichVolcanodata(upProteinsAD, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_AD/combinedData_AD_", xVars[i], "_upRegulated.csv"))
  downEnriched <- enrichVolcanodata(downProteinsAD, "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/combinedData_AD/combinedData_AD_", xVars[i], "_downRegulated.csv"))
}





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
knight_patients_CDR <- knight_plasma_patients[knight_plasma_patients$PA_DB_UID %in% knight_patients$PA_DB_UID, ] %>%
  select(CDR_at_blood_draw)
knight_patients_CDR <- knight_patients %>% mutate(CDR = knight_patients_CDR$CDR_at_blood_draw) %>%
  mutate(Age_CDR = avg_drawage * CDR)
knight_ratio <- knight_CSF_prots[-c(1:3)] / knight_plasma_prots[-c(1:3)]
knight_all_prots <- names(knight_ratio)
knight_interaction_lmodels <- generate_lmodels(knight_ratio, knight_patients_CDR, c("Sex", "Age_CDR"))
knight_interaction_lmsummary <- generate_lmsummary(knight_interaction_lmodels, knight_all_prots, c("Male", "Age*CDR"))
knight_interaction_volcanodata <- generate_volcanodata(knight_interaction_lmsummary, c("Male", "Age*CDR"))
generate_volcanoplot(knight_interaction_volcanodata, c("Male", "Age*CDR"), 2)

xVars <- c("Male", "Age*CDR")
for(i in 1:2) {
  upEnriched <- enrichVolcanodata(
    getTopProteins(knight_interaction_volcanodata[[i]], "up"), "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/Interaction_model/interaction_", xVars[i], "_upRegulated.csv"))
  upEnriched <- enrichVolcanodata(
    getTopProteins(knight_interaction_volcanodata[[i]], "down"), "ChEA_2022")$ChEA_2022
  write.csv(upEnriched, paste0(getwd(), "/Data/DiffAnalysis/Interaction_model/interaction_", xVars[i], "_downRegulated.csv"))
}


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

# Clustering and enrichment analysis 
ratio_AD_clustable <- table(ratio_AD_z_clustnum$cluster)
ratio_CO_clustable <- table(ratio_CO_z_clustnum$cluster)

ratio_AD_clust_enriched <- enrichClusts(ratio_AD_z_clustnum)
ratio_CO_clust_enriched <- enrichClusts(ratio_CO_z_clustnum)

for(i in 1:nclust) {
  write.csv(ratio_AD_clust_enriched[[i]], paste0(getwd(), "/Data/Cluster/combinedData_AD/cluster_", i, "_pathways.csv"))
  write.csv(ratio_CO_clust_enriched[[i]], paste0(getwd(), "/Data/Cluster/combinedData_CO/cluster_", i, "_pathways.csv"))
}

# differential sliding window analysis
sliding_window <- 10
# Key pathways across different windows --> i mean we can j do enrichR and see
# number of significant proteins --> idk how the natuer paper did it, cuz there
# isnt enough data for us to use a smaller date window
# -log10(q) of each protein --> new trajectory and cluster by varying significance


# combined Stanford volcano plot
# what proteins are correlating the most with albumin. 
# if the cluster score (cluster eigenvector? look into it.)
# include storage days for stanford by subtracting from today and rerun diff analysis
# introduce covariate for study design (binary variable for whether stanford/knight), add plasma_storagedays + drawdate_diff 
# sliding window 

write.csv(albumin_AD_enriched, "/home/sean777/albumin_AD_enriched.csv")




