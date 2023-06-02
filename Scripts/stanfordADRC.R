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

load_lib(c("tidyverse", "caret", "data.table", "igraph", "cluster", "purrr", 
           "corrplot", "dplyr", "ggplot2", "ggraph", "circlize", "Cairo", "ComplexHeatmap", 
           "ggrepel", "openxlsx", "stringr", "tidyr", "WGCNA"))
load_lib(c("corrr", "gprofiler2", "lsa", "factoextra", "pls", "loessclust"))

get_biodata <- function(path) {
  dir <- "/labs/twc/jarod/Data"
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
CSF_prots_fil_z <- scale(CSF_prots_fil[-c(1:3)])
plasma_prots_fil_z <- scale(plasma_prots_fil[-c(1:3)])
ratio <- CSF_prots_fil[-c(1:3)] / plasma_prots_fil[-c(1:3)]

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

# split dataset 
age_AD <- age_fil[patientdata_AD_index, ]
age_CO <- age_fil[patientdata_CO_index, ]
ratio_AD <- ratio[patientdata_AD_index, ]
ratio_CO <- ratio[patientdata_CO_index, ]
patientdata_fil <- patientdata_fil[patientdata_index, ]
patientdata_AD <- patientdata_fil[patientdata_fil$final_status == "AD", ]
patientdata_CO <- patientdata_fil[patientdata_fil$final_status == "CO", ]
ratio <- ratio[patientdata_index, ]
albumin_ratio <- ratio[, names(ratio) == "ALB.18380.78.3"]
ratio_z <- scale(ratio)
ratio_AD_z <- scale(ratio_AD)
ratio_CO_z <- scale(ratio_CO)
albumin_ratio_z <- scale(albumin_ratio)

# clustering for CO and DA
ratio_CO_z_dist <- dist(t(ratio_CO_z), method = "euclidean")
ratio_CO_z_clust <- hclust(ratio_CO_z_dist, method = "ward.D")
ratio_AD_z_dist <- dist(t(ratio_AD_z), method = "euclidean")
ratio_AD_z_clust <- hclust(ratio_AD_z_dist, method = "ward.D")

determine_numclust <- function(clust, dist) {
  silhouette <- vector()
  for (k in 2:10) {
    cut <- cutree(clust, k)
    vals <- silhouette(cut, dist)
    silhouette[k - 1] <- mean(vals[, 3])
  }
  plot(silhouette)
}

determine_numclust(ratio_CO_z_clust, ratio_CO_z_dist)
determine_numclust(ratio_AD_z_clust, ratio_AD_z_dist)
# we conclude both 10 clusters 

nclust <- 10

proteinByClust <- function(clust, dist, nclust) {
  clust <- cutree(clust, k = nclust) %>% 
    data.frame(cluster = .)
  clust["protein"] = rownames(clust)
  clust <- clust %>% relocate(cluster, .after = "protein")
  rownames(clust) <- c()
  return(clust)
}

ratio_AD_z_clustnum <- proteinByClust(ratio_AD_z_clust, ratio_AD_z_dist, nclust = 10)
ratio_CO_z_clustnum <- proteinByClust(ratio_CO_z_clust, ratio_CO_z_dist, nclust = 10)

proteinClustData <- function(clust, patient) {
  long <- cbind(data.frame(age = patient %>% dplyr::select(avg_drawage), 
                           ))
}

age = intake_age_seq), intake_pred_AD_z)
intake_pred_AD_z_long <- intake_pred_AD_z_long %>% gather(key = protein, value = value, -age) %>% 
  inner_join(intake_AD_clust, by = "protein") %>% 
  arrange(cluster, value, age)

generate_clusterplot <- function(data) {
  data %>% 
    ggplot(aes(x = age, y = value, group = protein, color = factor(cluster))) +
    geom_line(stat = "smooth", method = "loess", se = FALSE, linewidth = 0.5, alpha = 0.1) +
    scale_color_manual(
      values = c("#FF6F61", "#6B5B95", "#88B04B", "#FFA500", "#92A8D1", "#FF69B4", 
                 "#955251", "#008080", "#DA70D6", "#6A5ACD")) +
    labs(x = "X-axis", y = "Y-axis", color = "Cluster") +
    facet_wrap(~ cluster, ncol = 3, nrow = 4, labeller = labeller(cluster = as.character)) +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.8, alpha = 1, linetype = "solid", color = "white") +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.5, alpha = 0.8)
}




# LOOK AT RATIO BW PLASMA AND CSF ALBUMIN (albumin index). look at which ones 
# are correlated with albumin. (ALB.18380.78.3) Find which cluster it is in, etc 

# unified analysis: interaction model bw age and alzheimers disease status (extract CD? score)
# stratified analysis: remake the volcano plots 


# knight before p-adjusted VS stanford volcano plots --> compare enrichR pathways 

#interaction model bw age and alzheimers disease status (extract CD? score)



