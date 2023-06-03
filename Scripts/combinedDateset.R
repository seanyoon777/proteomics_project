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

match_by <- function(data, refdata) {
  data %>% 
    filter(PA_DB_UID %in% refdata$PA_DB_UID) %>% 
    arrange(match(PA_DB_UID, refdata$PA_DB_UID)) %>% 
    slice(order(match(PA_DB_UID, patientmeta_temp$PA_DB_UID)))
}

plasma_meta_temp <- match_by(plasma_meta_temp, plasma_patient_temp)
plasma_prots_temp <- match_by(plasma_prots_temp, plasma_patient_temp)
CSF_meta_temp <- match_by(CSF_meta_temp, CSF_patient_temp)
CSF_prots_temp <- match_by(CSF_prots_temp, CSF_patient_temp)

common_ID <- intersect(CSF_patient_temp$PA_DB_UID, plasma_patient_temp$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .)

CSF_meta_temp_fil <- match_by(CSF_meta_temp, common_ID)
CSF_prots_temp_fil <- match_by(CSF_prots_temp, common_ID)
CSF_patient_temp_fil <- match_by(CSF_patient_temp, common_ID)
plasma_meta_temp_fil <- match_by(plasma_meta_temp, common_ID)
plasma_prots_temp_fil <- match_by(plasma_prots_temp, common_ID)
plasma_patient_temp_fil <- match_by(plasma_patient_temp, common_ID)
patient_meta_temp_fil <- match_by(patientmeta_temp, common_ID)

patientdata_temp_fil <- merge(CSF_patient_temp_fil, plasma_patient_temp_fil, by = "PA_DB_UID") %>% 
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

date_window <- 1500
date_index <- patientdata_temp_fil$drawdate_diff <= date_window
patientdata_temp_fil_100 <- patientdata_temp_fil[date_index, ]
CSF_prots_temp_fil <- CSF_prots_temp_fil[date_index, ]
plasma_prots_temp_fil <- plasma_prots_temp_fil[date_index, ]

patientdata_temp <- patientdata_temp_fil_100 %>% 
  mutate(final_status = ifelse(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus), 
         age = abs(Plasma_drawage + CSF_drawage) / 2) %>% 
  select(PA_DB_UID, Sex, age, drawdate_diff, Plasma_drawdate, CSF_drawdate, final_status, storage_days) 


# split dataset 
patientdata_fil <- patientdata_fil[patientdata_index, ]
age_AD <- age_fil[patientdata_AD_index, ]
age_CO <- age_fil[patientdata_CO_index, ]
ratio <- ratio[patientdata_index, ]
ratio_AD <- ratio[patientdata_AD_index, ]
ratio_CO <- ratio[patientdata_CO_index, ]

AD_index_temp <- patientdata_temp$final_status == "AD"
CO_index_temp <- patientdata_temp$final_status == "CO"
ADCO_index <- AD_index_temp | CO_index_temp

CSF_prots_fil <- CSF_prots_fil[patientdata_index, ]
plasma_prots_fil <- plasma_prots_fil[patientdata_index, ]

all_index <- c(patientdata_index, ADCO_index)
all_AD_index <- c(patientdata_AD_index, AD_index_temp)
all_CO_index <- c(patientdata_CO_index, CO_index_temp)

all_patientdata <- bind_rows(patientdata_fil %>% mutate(Sex = Gender, drawdate_diff = abs(drawdate_diff)) %>% 
            dplyr::select(Sex, avg_drawage, final_status, drawdate_diff), 
          patientdata_temp %>% mutate(avg_drawage = age) %>% 
            dplyr::select(Sex, avg_drawage, final_status, drawdate_diff))
all_CSF <- bind_rows(CSF_prots_fil[-c(1:3)], CSF_prots_temp_fil[all_prots])
all_plasma <- bind_rows(plasma_prots_fil[-c(1:3)], plasma_prots_temp_fil[all_prots])

AD_index <- all_patientdata$final_status == "AD"
CO_index <- all_patientdata$final_status == "CO"
index <- AD_index | CO_index 

all_patientdata_AD <- all_patientdata[AD_index, ]
all_patientdata_CO <- all_patientdata[CO_index, ]
all_patientdata <- all_patientdata[index, ]
all_CSF <- all_CSF[index, ]
all_plasma <- all_plasma[index, ]

all_ratio <- all_CSF / all_plasma
all_ratio_AD <- all_ratio[all_patientdata$final_status == "AD", ]
all_ratio_CO <- all_ratio[all_patientdata$final_status == "CO", ]

all_ratio_z <- scale(all_ratio)
all_ratio_AD_z <- scale(all_ratio_AD)
all_ratio_CO_z <- scale(all_ratio_CO)

albumin_all_ratio <- all_ratio[, names(all_ratio) == "ALB.18380.78.3"]
albumin_all_ratio_z <- scale(albumin_all_ratio)


# Create histogram of population 
drawdate_windows <- 50 * (0:6)
drawdate_hists <- list()
for (i in 1:6) {
  patientdata_fil_bin <- all_patientdata[
    abs(all_patientdata$drawdate_diff) < drawdate_windows[i + 1] & 
      abs(all_patientdata$drawdate_diff) >= drawdate_windows[i], ]
  
  drawdate_hists[[i]] <- patientdata_fil_bin %>% 
    ggplot(aes(x = avg_drawage, fill = final_status)) + 
    geom_histogram(binwidth = 1, breaks = seq(30, 100, by = 10), 
                   position = "stack", show.legend = FALSE, color = "black") + 
    labs(title = paste0("Drawn within ", drawdate_windows[i], " and ", 
                        drawdate_windows[i+1], " days"), 
         x = "Average age when drawn", y = "Frequency") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
          plot.background = element_rect(fill = "white"),
          axis.title = element_text(size = 9.5),
          axis.text = element_text(size = 7),
          plot.title = element_text(size = 10, hjust = 0.5), 
          axis.line = element_blank()) +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    scale_fill_manual(values = c("CO" = "#00BFC4", "AD" = "#F8766D"))
  
  maxcount <- max(hist(patientdata_fil_bin$avg_drawage, 
                       breaks = seq(30, 100, by = 10))$counts)
  drawdate_hists[[i]] <- drawdate_hists[[i]] + 
    annotate("text", x = 95, y = maxcount, label = paste0("n = ", nrow(patientdata_fil_bin)), 
             size = 3)
}
#AD = red, CO = blue
gridExtra::grid.arrange(grobs = drawdate_hists, ncol = 3)



# predict data and trends 
age_min <- max(min(all_patientdata_AD$avg_drawage), min(all_patientdata_CO$avg_drawage))
age_max <- min(max(all_patientdata_AD$avg_drawage), max(all_patientdata_CO$avg_drawage))
age_seq <- seq(age_min, age_max, by = 0.25)

loess_predict <- function(x, Y, x_seq, xlab) {
  models <- vector("list", ncol(Y))
  data_pred <- data.frame(xlab = x_seq)
  for (i in 1:ncol(Y)) {
    models[[i]] <- loess(Y[, i] ~ x)
    data_pred[, i + 1] <- predict(models[[i]], x_seq)
  }
  colnames(data_pred) <- c(xlab, colnames(Y))
  return(data_pred)
}

predratio_AD_z <- loess_predict(all_patientdata_AD$avg_drawage, all_ratio_AD_z, age_seq, "Age")
predratio_CO_z <- loess_predict(all_patientdata_CO$avg_drawage, all_ratio_CO_z, age_seq, "Age")

# clustering for CO and DA
ratio_CO_z_dist <- dist(t(predratio_CO_z[-1]), method = "euclidean")
ratio_CO_z_clust <- hclust(ratio_CO_z_dist, method = "ward.D")
ratio_AD_z_dist <- dist(t(predratio_AD_z[-1]), method = "euclidean")
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
# we conclude both 8 clusters 

nclust <- 9

proteinByClust <- function(clust, dist, nclust) {
  clust <- cutree(clust, k = nclust) %>% 
    data.frame(cluster = .)
  clust["protein"] = rownames(clust)
  clust <- clust %>% relocate(cluster, .after = "protein")
  rownames(clust) <- c()
  return(clust)
}

ratio_AD_z_clustnum <- proteinByClust(ratio_AD_z_clust, ratio_AD_z_dist, nclust)
ratio_CO_z_clustnum <- proteinByClust(ratio_CO_z_clust, ratio_CO_z_dist, nclust) # tbh 8 clusters are better

proteinClustData <- function(data, clust) { 
  long <- data %>% gather(key = protein, value = value, -Age) %>% 
    inner_join(clust, by = "protein") %>% 
    arrange(cluster, Age, value)
  return(long)
}

ratio_AD_z_long <- proteinClustData(predratio_AD_z, ratio_AD_z_clustnum)
ratio_CO_z_long <- proteinClustData(predratio_CO_z, ratio_CO_z_clustnum)

generate_clusterplot <- function(data) {
  data %>% 
    ggplot(aes(x = Age, y = value, group = protein, color = factor(cluster))) +
    geom_line(stat = "smooth", method = "loess", se = FALSE, linewidth = 0.5, alpha = 0.1) +
    scale_color_manual(
      values = c("#FF6F61", "#6B5B95", "#88B04B", "#FFA500", "#92A8D1", "#FF69B4", 
                 "#955251", "#008080", "#DA70D6", "#6A5ACD")) +
    labs(x = "Age (years)", y = "CSF Plasma Protein Ratio (Z-Scored)", color = "Cluster") +
    facet_wrap(~ cluster, ncol = 3, nrow = 4, labeller = labeller(cluster = as.character)) +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.8, alpha = 1, linetype = "solid", color = "white") +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.5, alpha = 0.8) + 
    theme(strip.text = element_blank(), 
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_blank(),
          axis.line = element_blank(), 
          legend.position = "none") +
    labs(title = NULL)
}

generate_clusterplot(ratio_AD_z_long)
generate_clusterplot(ratio_CO_z_long)


# LOOK AT RATIO BW PLASMA AND CSF ALBUMIN (albumin index). look at which ones 
# are correlated with albumin. (ALB.18380.78.3) Find which cluster it is in, etc 
albumin_cluster_AD <- ratio_AD_z_clustnum[ratio_AD_z_clustnum$protein == "ALB.18380.78.3", ]$cluster
albumin_prots_AD <-ratio_AD_z_clustnum[ratio_AD_z_clustnum$cluster == albumin_cluster_AD, ]$protein
albumin_cluster_CO <- ratio_CO_z_clustnum[ratio_CO_z_clustnum$protein == "ALB.18380.78.3", ]$cluster
albumin_prots_CO <-ratio_CO_z_clustnum[ratio_CO_z_clustnum$cluster == albumin_cluster_CO, ]$protein
albumin_intersection <- intersect(albumin_prots_AD, albumin_prots_CO)

# unified analysis: interaction model bw age and alzheimers disease status (extract CD? score)


# stratified analysis: remake the volcano plots 
generate_lmodels <- function(protdata, patientdata) {
  lmodels <- list()
  for(i in 1:ncol(protdata)) {
    lmodels[[i]] <- summary(lm(protdata[, i] ~ patientdata$Sex + 
                                 patientdata$avg_drawage))
  }
  return(lmodels)
}

ratio_AD_models <- generate_lmodels(all_ratio_AD, all_patientdata_AD)
ratio_CO_models <- generate_lmodels(all_ratio_CO, all_patientdata_CO)

generate_lmsummary <- function(lmodels) {
  lmsummary <- data.frame()
  for (i in 1:length(all_prots)) {
    tidyresult <- data.frame(lmodels[[i]]$coefficients) %>% 
      select(Estimate, Std..Error, Pr...t..)
    colnames(tidyresult) <- c("Coefficient", "StdError", "Pval")
    tidyresult <- tidyresult[2:nrow(tidyresult), ] %>% 
      mutate(xVar = c("Male", "Age"), 
             Protein = all_prots[i]) %>% 
      dplyr::select(Protein, xVar, Coefficient, StdError, Pval)
    rownames(tidyresult) <- c()
    lmsummary <- rbind(lmsummary, tidyresult)
  }
  return(lmsummary)
}

ratio_AD_summary <- generate_lmsummary(ratio_AD_models)
ratio_CO_summary <- generate_lmsummary(ratio_CO_models)

xVars <- c("Male", "Age")
generate_volcanodata <- function(data) {
  volcano_data <- list()
  for (i in 1:length(xVars)) {
    volcano_data[[i]] <- data %>% filter(xVar == xVars[i]) %>% 
      mutate(padj = p.adjust(Pval, method = "fdr")) %>%
      mutate(qval = -log10(padj), 
             log2fc = Coefficient) %>% 
      select(Protein, Pval, padj, qval, log2fc) %>% 
      mutate(diffexpressed = ifelse(log2fc > 0 & qval >= 1,
                                    yes = "Upregulated", 
                                    no = ifelse(log2fc < 0 & qval >= 1,
                                                yes = "Downregulated", no="none")))
  }
  return(volcano_data)
}

ratio_AD_volcanodata <- generate_volcanodata(ratio_AD_summary)
ratio_CO_volcanodata <- generate_volcanodata(ratio_CO_summary)

generate_volcanodata_nonpadj <- function(data) {
  volcano_data <- list()
  for (i in 1:length(xVars)) {
    volcano_data[[i]] <- data %>% filter(xVar == xVars[i]) %>% 
      mutate(padj = Pval) %>%
      mutate(qval = -log10(padj), 
             log2fc = Coefficient) %>% 
      select(Protein, Pval, padj, qval, log2fc) %>% 
      mutate(diffexpressed = ifelse(log2fc > 0 & qval >= 1,
                                    yes = "Upregulated", 
                                    no = ifelse(log2fc < 0 & qval >= 1,
                                                yes = "Downregulated", no="none")))
  }
  return(volcano_data)
}

ratio_AD_volcanodata_nonpadj <- generate_volcanodata_nonpadj(ratio_AD_summary)

generate_volcanoplot <- function(volcanodata) {
  volcanoplot <- list()
  for(i in 1:length(xVars)) {
    top_genes <- head((volcanodata[[i]])[order(-volcanodata[[i]]$qval), ], 10)
    volcanoplot[[i]] <- ggplot(volcanodata[[i]], aes(x = log2fc, y = qval, color = factor(diffexpressed))) + 
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
      annotate("text", x = min((volcanodata[[i]])$log2fc)*1.1, y = max((volcanodata[[i]])$qval), 
               label = xVars[i], hjust = 0, size = 4) +
      geom_text_repel(data = top_genes, aes(label = Protein, vjust = qval, hjust = log2fc), 
                      size = 3, color = "black", min.segment.length = 0)
  }
  
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = 2)
}

generate_volcanoplot(ratio_AD_volcanodata)
generate_volcanoplot(ratio_CO_volcanodata)
generate_volcanoplot(ratio_AD_volcanodata_nonpadj)  


# Get enrichR pathways 
getTopProteins <- function(data, regulated = c("up", "down")) {
  list()
  for(i in 1:length(xVars)) {
    
  }
}

AD_top_prots <- 
  CO_top_prots <- age_top_up <- ratio_AD_volcanodata[[1]]
age_bottom_down 


# knight before p-adjusted VS stanford volcano plots --> compare enrichR pathways 

#interaction model bw age and alzheimers disease status (extract CD? score)



