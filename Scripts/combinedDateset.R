rm(list = ls())

source("Scripts/Helpers/init.R")
source("Scripts/Helpers/differentialAnalysis.R")
source("Scripts/Helpers/hierClust.R")
source("Scripts/Helpers/demoHist.R")
source("Scripts/Helpers/enrichAnalysis.R")
source("Scripts/Helpers/wgcna.R")

dir <- "~/Developer/Data/proteomics_project_data"
get_biodata <- function(path) {
  read.csv(paste(dir, path, sep = "/"))
}


# same thing with stanford data? ----
stanford_barcodes <- get_biodata("ADRC/barcodes.csv")
stanford_CSF_patients <- get_biodata("ADRC/CSF_patients.csv")
stanford_CSF_prots <- get_biodata("ADRC/CSF_prots.csv")
stanford_plasma_patients <- get_biodata("ADRC/plasma_patients.csv")
stanford_plasma_prots <- get_biodata("ADRC/plasma_prots.csv")

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
stanford_barcodes <- stanford_barcodes[stanford_gendermatch_index, ]

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
knight_plasma_meta <- get_biodata("Knight/plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
knight_plasma_prots <- get_biodata("Knight/plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
knight_plasma_patients <- read.xlsx(paste(dir, "Knight/plasma/Plasma_ BasicDemo.xlsx", sep = "/"))
knight_CSF_meta <- get_biodata("Knight/CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
knight_CSF_prots <- get_biodata("Knight/CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
knight_CSF_patients <- read.xlsx(paste(dir, "Knight/CSF/CSF_BasicDemo.xlsx", sep = "/"))

knight_patients_meta <- get_biodata("Knight/WashU_ADNI_demographic.csv")

knight_plasma_prots[4:ncol(knight_plasma_prots)] <- log10(knight_plasma_prots[4:ncol(knight_plasma_prots)])
knight_CSF_prots[4:ncol(knight_CSF_prots)] <- log10(knight_CSF_prots[4:ncol(knight_CSF_prots)])

protMeta <- read_csv("~/Developer/Data/proteomics_project_data/Knight/ProteinMetadata_with_LOD.csv", col_types = )
protMeta$SeqId = str_replace(protMeta$SeqId, "-", ".")
protMeta <- filter(protMeta, (Organism == "Human" | Organism == "HIV-1" | Organism == "HIV-2") & Type == "Protein")
all_prots_ID <- names(knight_plasma_prots)[4:ncol(knight_plasma_prots)] 
all_prots <- str_replace(all_prots_ID, ".*?(\\d+)\\.(\\d+)\\.\\d+$", "\\1.\\2")
my_order <- match(protMeta$SeqId, all_prots)
all_prots <- protMeta$Key_2
nprots <- length(all_prots)
my_order <- my_order + 3

knight_plasma_prots <- knight_plasma_prots[c(1:3, my_order)]
colnames(knight_plasma_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)

knight_CSF_prots <- knight_CSF_prots[c(1:3, my_order)]
colnames(knight_CSF_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)


all_prots <- intersect(names(knight_plasma_prots)[4:ncol(knight_plasma_prots)], 
                       names(knight_CSF_prots)[4:ncol(knight_CSF_prots)])
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
  mutate(DrawDate = convertToDateTime(as.numeric(DrawDate)))

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
all_plasma_AD <- all_plasma[all_patientdata$final_status == "AD", ]
all_plasma_CO <- all_plasma[all_patientdata$final_status == "CO", ]
all_CSF_AD <- all_CSF[all_patientdata$final_status == "AD", ]
all_CSF_CO <- all_CSF[all_patientdata$final_status == "CO", ]

all_ratio_z <- scale(all_ratio)
all_ratio_AD_z <- scale(all_ratio_AD)
all_ratio_CO_z <- scale(all_ratio_CO)
all_plasma_AD_z <- scale(all_plasma_AD)
all_plasma_CO_z <- scale(all_plasma_CO)
all_CSF_AD_z <- scale(all_CSF_AD)
all_CSF_CO_z <- scale(all_CSF_CO)

ALBUMIN_CODE <- "ALB.18380.78.3"
albumin_all_ratio <- all_ratio[, names(all_ratio) == ALBUMIN_CODE]
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
predplasma_AD_z <- loess_predict(all_patientdata_AD$avg_drawage, all_plasma_AD_z, age_seq, "Age")
predplasma_CO_z <- loess_predict(all_patientdata_CO$avg_drawage, all_plasma_CO_z, age_seq, "Age")
predCSF_AD_z <- loess_predict(all_patientdata_AD$avg_drawage, all_CSF_AD_z, age_seq, "Age")
predCSF_CO_z <- loess_predict(all_patientdata_CO$avg_drawage, all_CSF_CO_z, age_seq, "Age")

# clustering for CO and DA
ratio_CO_z_dist <- dist(t(predratio_CO_z[-1]), method = "euclidean")
ratio_CO_z_clust <- hclust(ratio_CO_z_dist, method = "complete")
ratio_AD_z_dist <- dist(t(predratio_AD_z[-1]), method = "euclidean")
ratio_AD_z_clust <- hclust(ratio_AD_z_dist, method = "complete")

plasma_CO_z_dist <- dist(t(predplasma_CO_z[-1]), method = "euclidean")
plasma_CO_z_clust <- hclust(plasma_CO_z_dist, method = "complete")
plasma_AD_z_dist <- dist(t(predplasma_AD_z[-1]), method = "euclidean")
plasma_AD_z_clust <- hclust(plasma_AD_z_dist, method = "complete")

CSF_CO_z_dist <- dist(t(predCSF_CO_z[-1]), method = "euclidean")
CSF_CO_z_clust <- hclust(CSF_CO_z_dist, method = "complete")
CSF_AD_z_dist <- dist(t(predCSF_AD_z[-1]), method = "euclidean")
CSF_AD_z_clust <- hclust(CSF_AD_z_dist, method = "complete")

determine_numclust(ratio_CO_z_clust, ratio_CO_z_dist)
determine_numclust(ratio_AD_z_clust, ratio_AD_z_dist)

determine_numclust(CSF_CO_z_clust, CSF_CO_z_dist)
determine_numclust(CSF_AD_z_clust, CSF_AD_z_dist)

determine_numclust(plasma_CO_z_clust, plasma_CO_z_dist)
determine_numclust(plasma_AD_z_clust, plasma_AD_z_dist)
# we conclude both 9 clusters 

nclust <- 10

ratio_AD_z_clustnum <- proteinByClust(ratio_AD_z_clust, ratio_AD_z_dist, nclust)
ratio_CO_z_clustnum <- proteinByClust(ratio_CO_z_clust, ratio_CO_z_dist, nclust) # tbh 8 clusters are better

ratio_AD_z_long <- proteinClustData(predratio_AD_z, ratio_AD_z_clustnum)
ratio_CO_z_long <- proteinClustData(predratio_CO_z, ratio_CO_z_clustnum)

generate_clusterplot(ratio_AD_z_long)
generate_clusterplot(ratio_CO_z_long)
ratio_CO_z_clustnum[ratio_CO_z_clustnum$cluster == 10, ]

nclust <- 9
CSF_AD_z_clustnum <- proteinByClust(CSF_AD_z_clust, CSF_AD_z_dist, nclust)
CSF_CO_z_clustnum <- proteinByClust(CSF_CO_z_clust, CSF_CO_z_dist, nclust)
CSF_AD_z_long <- proteinClustData(predCSF_AD_z, CSF_AD_z_clustnum)
CSF_CO_z_long <- proteinClustData(predCSF_CO_z, CSF_CO_z_clustnum)
generate_clusterplot(CSF_AD_z_long)
generate_clusterplot(CSF_CO_z_long)
CSF_AD_z_clustnum[CSF_AD_z_clustnum$cluster == 9,]

nclust <- 8
plasma_AD_z_clustnum <- proteinByClust(plasma_AD_z_clust, plasma_AD_z_dist, nclust)
plasma_CO_z_clustnum <- proteinByClust(plasma_CO_z_clust, plasma_CO_z_dist, nclust)
plasma_AD_z_long <- proteinClustData(predplasma_AD_z, plasma_AD_z_clustnum)
plasma_CO_z_long <- proteinClustData(predplasma_CO_z, plasma_CO_z_clustnum)
generate_clusterplot(plasma_AD_z_long)
generate_clusterplot(plasma_CO_z_long)
plasma_CO_z_clustnum[plasma_CO_z_clustnum$cluster == 7,]

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
generate_volcanoplot(ratio_AD_volcanodata, xLabs, 2, 10, c(FALSE, FALSE))
generate_volcanoplot(ratio_CO_volcanodata, xLabs, 2, 10, c(FALSE, FALSE))
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
generate_volcanoplot(ratio_volcanodata, c("Male", "Age", "AD"), 3, 10, c(FALSE, FALSE, FALSE))

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

# Just curious how these are different. 
all_CSF_lmodels <- generate_lmodels(all_CSF, all_patientdata, 
  c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect"))
all_CSF_lmsummary <- generate_lmsummary(all_CSF_lmodels, all_prots, c("Male", "Age", "AD", "drawdate_diff", "storage_days", "batch_effect"))
all_CSF_volcanodata <- generate_volcanodata(all_CSF_lmsummary, c("Male", "Age", "AD"))
generate_volcanoplot(all_CSF_volcanodata, c("Male", "Age", "AD"), 3, 10, c(FALSE, FALSE, TRUE))

all_plasma_lmodels <- generate_lmodels(all_plasma, all_patientdata, 
                                    c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect"))
all_plasma_lmsummary <- generate_lmsummary(all_plasma_lmodels, all_prots, c("Male", "Age", "AD", "drawdate_diff", "storage_days", "batch_effect"))
all_plasma_volcanodata <- generate_volcanodata(all_plasma_lmsummary, c("Male", "Age", "AD"))
generate_volcanoplot(all_plasma_volcanodata, c("Male", "Age", "AD"), 3, 10, c(FALSE, FALSE, TRUE))

# knight before p-adjusted VS stanford volcano plots --> compare enrichR pathways 
knight_ratio <- knight_CSF_prots[-c(1:3)] / knight_plasma_prots[-c(1:3)]
knight_all_prots <- colnames(knight_ratio)

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

generate_volcanoplot(knight_AD_volcanodata, c("Male", "Age"), 2, 10, c(FALSE, FALSE))
generate_volcanoplot(knight_CO_volcanodata, c("Male", "Age"), 2, 10, c(FALSE, FALSE))

#interaction model bw age and alzheimers disease status (extract CD? score)
knight_patients_CDR <- knight_plasma_patients[knight_plasma_patients$PA_DB_UID %in% knight_patients$PA_DB_UID, ] %>%
  select(CDR_at_blood_draw)
knight_patients_CDR <- knight_patients %>% mutate(CDR = knight_patients_CDR$CDR_at_blood_draw)
knight_patients_CDR$CDR <- if_else(is.na(knight_patients_CDR$CDR) & knight_patients_CDR$final_status == "CO", 
                                   0.00, knight_patients_CDR$CDR)
knight_patients_CDR <- knight_patients_CDR %>% mutate(Age_CDR = avg_drawage * CDR)
knight_ratio <- knight_CSF_prots[-c(1:3)] / knight_plasma_prots[-c(1:3)]
knight_all_prots <- names(knight_ratio)
knight_interaction_lmodels <- generate_lmodels(knight_ratio, knight_patients_CDR, c("Sex", "Age_CDR", "drawdate_diff", "storage_days"))
knight_interaction_lmsummary <- generate_lmsummary(knight_interaction_lmodels, knight_all_prots, c("Male", "Age*CDR", "drawdate_diff", "storage_days"))
knight_interaction_volcanodata <- generate_volcanodata(knight_interaction_lmsummary, c("Male", "Age*CDR", "drawdate_diff", "storage_days"))
generate_volcanoplot(knight_interaction_volcanodata, c("Male", "Age*CDR"), 2, 10, c(FALSE, FALSE))

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


# what proteins are correlating the most with albumin. 
# if the cluster score (cluster eigenvector? look into it.)
# sliding window 

euclidean <- function(a, b) sqrt(sum((a - b)^2))
cosine <- function(a, b) (sum(a * b)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
predratio_AD_znorm <- apply(predratio_AD_z[-1], 2, function(x) (x - min(x)) / (max(x) - min(x)))
predratio_CO_znorm <- apply(predratio_CO_z[-1], 2, function(x) (x - min(x)) / (max(x) - min(x)))
all_prots <- names(predratio_CO_z[-1])
ratio_similarities <- data.frame(protein = all_prots)
ratio_similarities <- inner_join(data.frame(protein = all_prots), ratio_AD_z_clustnum)
ratio_similarities <- inner_join(ratio_similarities, ratio_CO_z_clustnum, by = join_by(protein))
albumin_similarities <- ratio_similarities
nprots <- length(all_prots)
for (i in 1:nprots) {
  ratio_similarities[i, 4] <- euclidean(predratio_AD_z[, i + 1], predratio_CO_z[, i + 1]) 
  ratio_similarities[i, 5] <- cosine(predratio_AD_z[, i + 1], predratio_CO_z[, i + 1]) 
}
colnames(ratio_similarities) <- c("protein", "AD_cluster", "CO_cluster", "euclidean_dist", "cosine_sim")
for (i in 1:nprots) {
  ratio_similarities[i, 6] <- euclidean(predratio_AD_znorm[, i], predratio_CO_znorm[, i]) 
  ratio_similarities[i, 7] <- cosine(predratio_AD_znorm[, i], predratio_CO_znorm[, i])  
}
colnames(ratio_similarities) <- c("protein", "AD_cluster", "CO_cluster", "euclidean", "cosine", 
                                  "norm_euclidean", "norm_cosine")

write.csv(ratio_similarities, paste0(getwd(), "/Data/Trajectories/ratioSimilarities.csv"))

predratio_AD_znorm <- data.frame(predratio_AD_znorm)
predratio_CO_znorm <- data.frame(predratio_CO_znorm)
for(i in 1:(nprots - 1)) {
  albumin_similarities[i, 4] <- euclidean(predratio_AD_z[, i + 1], predratio_AD_z[ALBUMIN_CODE]) 
  albumin_similarities[i, 5] <- cosine(predratio_AD_z[, i + 1], predratio_AD_z[ALBUMIN_CODE]) 
  albumin_similarities[i, 6] <- euclidean(predratio_AD_znorm[, i], predratio_AD_znorm[ALBUMIN_CODE]) 
  albumin_similarities[i, 7] <- cosine(predratio_AD_znorm[, i], predratio_AD_znorm[ALBUMIN_CODE])  
  albumin_similarities[i, 8] <- euclidean(predratio_CO_z[, i + 1], predratio_CO_z[ALBUMIN_CODE]) 
  albumin_similarities[i, 9] <- cosine(predratio_CO_z[, i + 1], predratio_CO_z[ALBUMIN_CODE]) 
  albumin_similarities[i, 10] <- euclidean(predratio_CO_znorm[, i], predratio_CO_znorm[ALBUMIN_CODE]) 
  albumin_similarities[i, 11] <- cosine(predratio_CO_znorm[, i], predratio_CO_znorm[ALBUMIN_CODE])  
}

colnames(albumin_similarities) <- c("protein", "AD_cluster", "CO_cluster", "AD_euclidean", 
                                    "AD_cosine", "AD_norm_euclidean", "AD_norm_cosine", 
                                    "CO_euclidean", "CO_cosine", "CO_norm_euclidean", 
                                    "CO_norm_cosine")

albumin_similarities <- albumin_similarities[albumin_similarities$protein != ALBUMIN_CODE, ]
write.csv(albumin_similarities, paste0(getwd(), "/Data/Trajectories/albuminSimilarities.csv"))

# Ratio stability wrt blood/CSF drawdate
## 1. DEA for each bin
ratio_stability_DEA <- function(all_patientdata, all_ratio, interval) {
  num <- 120 / interval
  drawdate_bins <- c(0, interval * 0:num)
  volcanodata <- vector(mode = 'list', length = num)
  stability_volcanoplot <- vector(mode = 'list', length=num)
  
  for(i in 1:(num + 1)) {
    idx <- all_patientdata$drawdate_diff <= drawdate_bins[i + 1] & all_patientdata$drawdate_diff >= drawdate_bins[i]
    new_patientdata <- all_patientdata[idx, ]
    new_ratiodata <- all_ratio[idx, ]
    
    if(sum(new_patientdata$batch_effect == 1) * sum(new_patientdata$batch_effect == 0) == 0) {
      lmodels <- generate_lmodels(new_ratiodata, new_patientdata, c("Sex", "avg_drawage", "final_status"))
      lmsummary <- generate_lmsummary(lmodels, all_prots, c("Male", "Age", "AD"))
    } else {
      lmodels <- generate_lmodels(new_ratiodata, new_patientdata, c("Sex", "avg_drawage", "final_status", "batch_effect"))
      lmsummary <- generate_lmsummary(lmodels, all_prots, c("Male", "Age", "AD", "batch_effect"))
    }
    volcanodata[[i]] <- generate_volcanodata(lmsummary, c("Male", "Age", "AD"))
  } 
  
  return(volcanodata)
}
ratio_stability_volcanodata_0 <- ratio_stability_DEA(all_patientdata, all_ratio, interval = 20)
ratio_stability_volcanodata_40_0 <- ratio_stability_DEA(all_patientdata, all_ratio, interval = 40)


ratio_stability_plot <- function(all_patientdata, ratio_stability_volcanodata, interval, num_top_proteins, boxpadding, factor_num, ncol = 4) {
  titles <- c('Sex', 'Aging', 'Alzheimers')
  volcanoplot <- list()
  num <- 120 / interval 
  drawdate_bins <- c(0, interval * 0:num+ 0.01)
  for(i in 1:(num + 1)) {
    nsamples <- sum(all_patientdata$drawdate_diff <= drawdate_bins[i + 1] & all_patientdata$drawdate_diff >= drawdate_bins[i])
    volcanodata_temp <- ratio_stability_volcanodata[[i]][[factor_num]]
    top_volcanodata <- volcanodata_temp[volcanodata_temp$diffexpressed != "none", ]
    top_genes <- head(top_volcanodata[order(-top_volcanodata$qval), ], num_top_proteins)
    max_qval <- max(abs(volcanodata_temp$qval), na.rm = TRUE)
    
    volcanodata_temp$point_size <- abs(volcanodata_temp$qval)^2 / max_qval
    if (sum(volcanodata_temp$diffexpressed != "none") == 0) {
      volcanodata_temp$point_size <- 0.1
      scale_size <- c(0.5, 1.5)
    } else {
      scale_size <- c(1, 6)
    }
    
    volcanoplot[[i]] <- ggplot(volcanodata_temp, aes(x = log2fc, y = qval, color = factor(diffexpressed))) + 
      geom_point(aes(size = point_size, alpha = abs(qval)/max_qval), na.rm = T) +
      scale_size_continuous(range = scale_size) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      #geom_point(color = "black", shape = 21, size = 0.5, stroke = 0.5, na.rm = T) +
      #geom_text_repel(max.overlaps = 10, aes(label = delabel)) + 
      theme_bw(base_size = 16) +
      theme(legend.position = "none") +
      #ggtitle(label = paste(str_split(plot_title, '_', simplify = T), sep = "", collapse = " ")) + 
      xlab("log2 FC") +
      ylab(expression(-log[10]("q"))) +
      scale_color_manual(values = c("Upregulated" = "indianred1",
                                    "Downregulated" = "royalblue1",
                                    "none" = "#2c2c2c")) +
      annotate("text", x = min((volcanodata_temp)$log2fc)*1.1, y = max((volcanodata_temp)$qval), 
               label = paste0(toString(floor(drawdate_bins[i])), '-',toString(floor(drawdate_bins[i + 1])), 
                              ' (', nsamples,' samples)'),
               hjust = 0, size = 4) +
      geom_text_repel(data = top_genes, aes(label = Protein, y = qval, x = log2fc),
                      size = 3, color = "black", nudge_y = -0.2) + 
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
            plot.background = element_rect(fill = "white"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            plot.title = element_text(size = 10, hjust = 0.5), 
            axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(family = "Arial")) 
  }
  
  title <- paste0(titles[factor_num], " Proteome by Drawdate Difference")
  title_grob <- textGrob(title, gp=gpar(fontsize=20, fontface="bold"))
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = ncol, top = title_grob)
}
ratio_stability_plot(all_patientdata, ratio_stability_volcanodata_0, num_top_proteins = 12, 20, boxpadding = .5, factor_num = 3)
ratio_stability_plot(all_patientdata, ratio_stability_volcanodata_40_0, num_top_proteins = 12, 40, boxpadding = .5, factor_num = 1, ncol = 2)

padj_rank <- function() {
  
}
#Use p values to plot the proteins. Rank P value statistic for each protein. 









#
geom_text_repel(data = top_genes, aes(label = Protein, vjust = qval, hjust = log2fc),
                size = 3, color = "black", box.padding = boxpadding, 
                max.overlaps = Inf, force_pull = 20, max.iter = 50000) 




inverse_normal_transform <- function(x) {
  ranks <- rank(x)
  percentiles <- ranks / (length(x) + 1)
  transformed <- qnorm(percentiles)
  return(transformed)
}
all_ratio_int <- data.frame(apply(all_ratio, 2, inverse_normal_transform))
ratio_int_stability_volcanodata <- ratio_stability_DEA(all_patientdata, all_ratio_int, interval = 40)
ratio_stability_plot(all_patientdata, ratio_int_stability_volcanodata, num_top_proteins = 12, 40, boxpadding = .5, factor_num = 3)



#suddenly thought of a scary idea
stanford_ratio <- stanford_CSF_prots[-c(1:3)] / stanford_plasma_prots[-c(1:3)]
stanford_all_prots <- colnames(stanford_ratio)
stanford_lmodel <- generate_lmodels(stanford_ratio, stanford_patients, c("avg_drawage", "Gender", "final_status"))
stanford_lmsummary <- generate_lmsummary(stanford_lmodel, stanford_all_prots, c("Age", "Male", "AD"))
stanford_volcanodata <- generate_volcanodata(stanford_lmsummary, c("Age", "Male", "AD"))
generate_volcanoplot(stanford_volcanodata, c("Age", "Male", "AD"), 3, 12)

knight_ratio <- knight_CSF_prots[-c(1:3)] / knight_plasma_prots[-c(1:3)]
knight_all_prots <- colnames(knight_ratio)
knight_lmodel <- generate_lmodels(knight_ratio, knight_patients, c("avg_drawage", "Sex", "final_status"))
knight_lmsummary <- generate_lmsummary(knight_lmodel, knight_all_prots, c("Age", "Male", "AD"))
knight_volcanodata <- generate_volcanodata(knight_lmsummary, c("Age", "Male", "AD"))
generate_volcanoplot(knight_volcanodata, c("Age", "Male", "AD"), 3)
# almost fucked myself over


decide_soft_threshold(all_plasma) # threshold = 6
decide_soft_threshold(all_CSF)  # threshold = 6 - 8?
decide_soft_threshold(all_ratio)  # threshold >= 2, but this shoudl be true since ratio X have dims. 

plasma_net <- cluster_dendrogram(all_plasma, 6, "Plasma")
CSF_net <- cluster_dendrogram(all_CSF, 6, "CSF")
ratio_net <- cluster_dendrogram(all_ratio, 2, "Ratio")

cluster_heatmap_dendrogram(all_plasma, plasma_net, "Plasma", 6)
cluster_heatmap_dendrogram(all_CSF, CSF_net, "CSF", 6)
cluster_heatmap_dendrogram(all_ratio, ratio_net, "Ratio", 2)

module_membership <- data.frame(
  Gene = names(all_ratio), 
  ModuleColor = plasma_net$colors
)
write.csv(module_membership, file = "generated_data/plasma_module_membership.csv", row.names = FALSE)


heatmap_raw <- function() {
  
}

heatmap_raw_wgcna <- function() {
  
}
corr_matrix <- cor(all_ratio)
hc <- hclust(as.dist(1 - corr_matrix), method = "average")
order_of_genes <- order.dendrogram(as.dendrogram(hc))
corr_matrix2 <- corr_matrix[order_of_genes, order_of_genes]
melt_corr_matrix2 <- melt(corr_matrix2)

geneTree = ratio_net$dendrograms[[1]]

# Order of genes as they appear in the dendrogram
order_of_genes = order.dendrogram(as.dendrogram(geneTree))

corr_matrix <- corr_matrix[order_of_genes, order_of_genes]
corr_matrix <- melt(corr_matrix)

corr_matrix %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colors = colorRampPalette(c("red", "orange", "yellow", "white"))(200), guide = "none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())



corrplot(plasma_corr_matrix,
         col = colorRampPalette(c("red", "orange", "yellow", "white"))(200),
         addCoef.col = "black", cl.pos = 'n')

c("Male", "Age", "AD")



dendrogram <- hclust(dist_matrix)
cophenetic_matrix <- cophenetic(dendrogram)
correlation <- cor(as.vector(dist_matrix), as.vector(cophenetic_matrix))




Comparing two trees, particularly in the context of computational biology or computer science, can be done through various methods depending on the specific need and the type of trees. Here are some common methods:
  
  Tree Edit Distance: This method calculates how many operations (insertions, deletions, substitutions) are required to change one tree into the other. Algorithms like the Zhang-Shasha algorithm can be used to compute this.
Robinson-Foulds Metric: Often used to compare phylogenetic trees, the Robinson-Foulds metric counts the number of splits or clades that are found in one tree but not in the other.
Subtree Pruning and Regrafting (SPR) Distance: The SPR distance measures the minimum number of SPR moves required to transform one tree into the other. An SPR move involves pruning a subtree and regrafting it elsewhere in the tree.
Tanglegram: A tanglegram is a visual method for comparing two trees. The leaves of the trees, which must have the same set of labels, are drawn in two parallel lines, with lines connecting matching labels. Crossing lines indicate discrepancies between the trees.
Quartet Distance: This method looks at all the possible sets of four leaves and checks how often the two trees agree or disagree on the topology of these quartets.
Triplets Distance: Similar to the Quartet Distance but considers sets of three leaves.
Normalized Compression Distance (NCD): A method based on Kolmogorov complexity that can be used to define a similarity metric between two trees.
Jaccard Similarity: If the trees can be represented as sets (e.g., sets of edges or paths), the Jaccard similarity coefficient can be used to measure the similarity between these sets.





MEs <- moduleEigengenes(all_plasma, moduleColors)$eigengenes
weight = as.data.frame(all_patientdata$final_status);
weight <- as.data.frame(ifelse(all_patientdata$final_status == 'CO', 0, 1))

names(weight) = "weight"

MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(all_plasma);
nSamples = nrow(all_plasma);

dissTOM = 1-TOMsimilarityFromExpr(all_plasma, power = 6);
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


# Recalculate module eigengenes
MEs = moduleEigengenes(all_plasma, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes





# Recalculate MEs with color labels
MEs0 = moduleEigengenes(all_plasma, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits <- all_patientdata[c("Sex", "avg_drawage", "final_status")]
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);





ADJ_plasma <- abs(cor(all_plasma, use = "p"))^11
k = as.vector(apply(ADJ_plasma, 2, sum, na.rm = T))
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")



library(ggplot2)
library(reshape2)

# Suppose df is your data frame with multiple columns
# Melt the data to long format
all_ratio_melt <- melt(all_ratio)

new_age <- rep(all_patientdata$avg_drawage, times = ncol(all_ratio))

all_ratio_melt <- all_ratio_melt %>% mutate(age = new_age)

ggplot(data=all_ratio_melt,
       aes(x=age, y=value, colour=variable)) +
  geom_line(alpha = 0.1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "X-axis label", y = "Y-axis label", color = "Legend title")


