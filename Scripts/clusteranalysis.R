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

load_lib(c("tidyverse", "caret", "data.table", "corrr", "igraph", "cluster", "factoextra", "purrr", 
           "corrplot", "dplyr", "ggplot2", "pls", "ggraph", "circlize", "Cairo", "ComplexHeatmap", 
           "loessclust", "ggrepel", "lsa"))

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

intake_fil <- 1 - CSF_prot_fil[4:ncol(CSF_prot_fil)] / plasma_prot_fil[4:ncol(plasma_prot_fil)]

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

intake_strict_lm_models <- list()
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


eliminate_outliers <- function(data, var) {
  vec <- as.numeric(data[[var]])
  ub <- mean(vec) + IQR(vec) * 1.5
  lb <- mean(vec) - IQR(vec) * 1.5
  index <- vec <= ub | vec >= lb
  data[index, ]
}

generate_volcanodata <- function(data, variable) {
  volcano_data <- data %>% filter(xVar == variable) %>% 
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

generate_volcanodata_outliers <- function(data, variable) {
  generate_volcanodata(eliminate_outliers(data, "Coefficient"), variable)
}

generate_volcanoplot <- function(volcano_data, variable) {
  top_genes <- head(volcano_data[order(-volcano_data$qval), ], 10)
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
             label = variable, hjust = 0, size = 4) +
    geom_text_repel(data = top_genes, aes(label = Protein, vjust = qval, hjust = log2fc), 
                     size = 3, color = "black", min.segment.length = 0)
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

# write csv files 
for (i in 1:6) {
  write.csv(knight_strict_volcanodata[[i]], 
            paste0("data/generated_data/knight_volcanodata/knight_strict_volcanodata_", xVars[i]))
}

for (i in 1:6) {
  write.csv(knight_lenient_volcanodata[[i]], 
            paste0("data/generated_data/knight_volcanodata/knight_lenient_volcanodata_", xVars[i]))
}


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

for (i in 1:6) {
  write.csv(knight_lenient_volcanodata[[i]], 
            paste0("data/generated_data/knight_lenient_volcanodata_", xVars[i]))
}


# CLUSTERING ANALYSIS
AD_index <- patientdata_summary$final_status == "AD"
patientdata_AD <- patientdata_summary[AD_index, ]
CO_index <- patientdata_summary$final_status == "CO"
patientdata_CO <- patientdata_summary[CO_index, ]

intake_AD <- intake_fil[AD_index, ]
intake_CO <- intake_fil[CO_index, ]

age_min <- max(min(patientdata_AD$age), min(patientdata_CO$age))
age_max <- min(max(patientdata_AD$age), max(patientdata_CO$age))
intake_age_seq <- seq(age_min, age_max, by = 0.25)


generate_loess_preddata <- function(x, Y, x_seq, xlab) {
  models <- vector("list", ncol(Y))
  data_pred <- data.frame(xlab = x_seq)
  for (i in 1:ncol(Y)) {
    models[[i]] <- loess(Y[[i]] ~ x)
    data_pred[, i + 1] <- predict(models[[i]], x_seq)
  }
  colnames(data_pred) <- c(xlab, all_prots)
  return(data_pred)
}


intake_pred_AD <- generate_loess_preddata(patientdata_AD$age, intake_AD, intake_age_seq, "Age")
intake_pred_AD_z <- scale(intake_pred_AD[-1])
intake_pred_CO <- generate_loess_preddata(patientdata_CO$age, intake_CO, intake_age_seq, "Age")
intake_pred_CO_z <- scale(intake_pred_CO[-1])

intake_AD_z_dist <- dist(t(intake_pred_AD_z), method = "euclidean")
intake_AD_z_clust <- hclust(intake_AD_z_dist, method = "ward.D")
intake_CO_z_dist <- dist(t(intake_pred_CO_z), method = "euclidean")
intake_CO_z_clust <- hclust(intake_CO_z_dist, method = "ward.D")

silhouette_list <- vector()
for (k in 2:10) {
  AD_cut_clusters <- cutree(intake_AD_z_clust, k)
  AD_vals <- silhouette(AD_cut_clusters, intake_AD_z_dist)
  silhouette_val <- mean(AD_vals[, 3])
  silhouette_list[k - 1] <- silhouette_val
}
plot(silhouette_list)
# 6 or 9 clusters seems optimal? lets check for CO. 

silhouette_list_CO <- vector()
for (k in 2:10) {
  CO_cut_clusters <- cutree(intake_CO_z_clust, k)
  CO_vals <- silhouette(CO_cut_clusters, intake_CO_z_dist)
  silhouette_val <- mean(CO_vals[, 3])
  silhouette_list_CO[k - 1] <- silhouette_val
}
plot(silhouette_list_CO)
# 6 or 10 clusters seem optimal. 
# so lets try 6 for both? 
nclust <- 10

intake_AD_clust <- cutree(intake_AD_z_clust, k = nclust) %>% 
  data.frame(cluster = .)
intake_AD_clust["protein"] = rownames(intake_AD_clust)
intake_AD_clust <- intake_AD_clust %>% select(protein, cluster)
rownames(intake_AD_clust) <- c()

intake_CO_clust <- cutree(intake_CO_z_clust, k = nclust) %>%
  data.frame(cluster = .)
intake_CO_clust["protein"] = rownames(intake_CO_clust)
intake_CO_clust <- intake_CO_clust %>% select(protein, cluster)
rownames(intake_CO_clust) <- c()


# lets see which ones changed clusters 
sum(intake_AD_clust$cluster == intake_CO_clust$cluster) # theres a lot of proteins that changed clusters (6017)

# we need a different method...
# maybe first try plotting? 
intake_pred_AD_z_long <- cbind(data.frame(age = intake_age_seq), intake_pred_AD_z)
intake_pred_AD_z_long <- intake_pred_AD_z_long %>% gather(key = protein, value = value, -age) %>% 
  inner_join(intake_AD_clust, by = "protein") %>% 
  arrange(cluster, value, age)

intake_pred_CO_z_long <- cbind(data.frame(age = intake_age_seq), intake_pred_CO_z)
intake_pred_CO_z_long <- intake_pred_CO_z_long %>% gather(key = protein, value = value, -age) %>% 
  inner_join(intake_CO_clust, by = "protein") %>% 
  arrange(cluster, protein, age) 



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

#interesting plots...
generate_clusterplot(intake_pred_CO_z_long)
generate_clusterplot(intake_pred_AD_z_long)

# tried to see macro trend change, but not sure
# let's see specific proteins. choose a protein and plot AD and CO 
AD_CO_comparison <- function(protein) {
  AD_CO_data <- data.frame(age = intake_age_seq, 
                           CO = intake_pred_CO_z[, protein], AD = intake_pred_AD_z[, protein]) %>%
    gather(key = diagnosis, value = value, -age)
  ggplot(AD_CO_data, aes(age, value)) + 
    geom_line(aes(color = diagnosis))
}

AD_CO_comparison(all_prots[4])

# quantify how different each of the trajectories are and rank them?
# 1. maybe not on values, but how the "shape" changed? in this case, KS test
# 2. maybe on values 

euclidean <- function(a, b) sqrt(sum((a - b)^2))
intake_pred_AD_znorm <- apply(intake_pred_AD_z, 2, function(x) (x - min(x)) / (max(x) - min(x)))
intake_pred_CO_znorm <- apply(intake_pred_CO_z, 2, function(x) (x - min(x)) / (max(x) - min(x)))
intake_traj_dist <- data.frame(protein = all_prots)
intake_traj_dist <- inner_join(data.frame(protein = all_prots), intake_AD_clust)
intake_traj_dist <- inner_join(intake_traj_dist, intake_CO_clust, by = join_by(protein))
for (i in 1:nprots) {
  intake_traj_dist[i, 4] <- euclidean(intake_pred_AD_z[, i], intake_pred_CO_z[, i]) 
  intake_traj_dist[i, 5] <- cosine(intake_pred_AD_z[, i], intake_pred_CO_z[, i]) 
}
colnames(intake_traj_dist) <- c("protein", "AD_cluster", "CO_cluster", "euclidean_dist", "cosine_sim")
for (i in 1:nprots) {
  intake_traj_dist[i, 6] <- euclidean(intake_pred_AD_znorm[, i], intake_pred_CO_znorm[, i]) 
  intake_traj_dist[i, 7] <- cosine(intake_pred_AD_znorm[, i], intake_pred_CO_znorm[, i])  
}
colnames(intake_traj_dist) <- c("protein", "AD_cluster", "CO_cluster", "euclidean", "cosine", 
                                "norm_euclidean", "norm_cosine")

write.csv(intake_traj_dist, "data/generated_data/knight_intake_traj_dist.csv")




# same thing with stanford data? ----
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

patientdata_filtered <- plasma_patient2 %>% 
  rename(Plasma_Barcode = Barcode, Plasma_Zscore = ConnectivityZscore) %>% 
  filter(Plasma_Barcode %in% barcodes2$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes2$Plasma_Barcode)) 

stanford_gendermatch_index <- CSF_patient2_fil$Gender == patientdata_filtered$Gender
CSF_patient2_fil <- CSF_patient2_fil[stanford_gendermatch_index, ]
patientdata_filtered <- patientdata_filtered[stanford_gendermatch_index, ]
age2_filtered <- age2_filtered[stanford_gendermatch_index, ]
CSF_protein2_fil <- CSF_protein2_fil[stanford_gendermatch_index, ]
plasma_protein2_fil <- plasma_protein2_fil[stanford_gendermatch_index, ]

patientdata2_fil <- patientdata_filtered %>% 
  mutate(CSF_Barcode = CSF_patient2_fil$Barcode, Plasma_drawage = Age, CSF_drawage = CSF_patient2_fil$Age, 
         Plasma_storagedays = Storage_days, CSF_storagedays = CSF_patient2_fil$Storage_days,
         CSF_days = CSF_patient2_fil$Storage_days, Plasma_drawdate = strptime(Date.of.draw, format = "%m/%d/%y"), 
         CSF_drawdate = strptime(CSF_patient2_fil$Date_of_CSF_sample, format = "%m/%d/%y"), 
         Plasma_drawstatus = Diagnosis_group, CSF_drawstatus = CSF_patient2_fil$Diagnosis_group) 
patientdata2_fil <- patientdata2_fil %>% 
  group_by(Plasma_Barcode) %>%
  mutate(avg_drawage = (Plasma_drawage + as.numeric(CSF_drawage)) / 2, 
         drawdate_diff = difftime(Plasma_drawdate, CSF_drawdate, units = "days"), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>%
  select(Plasma_Barcode, CSF_Barcode, Gender, avg_drawage, Plasma_drawdate, Plasma_drawage, Plasma_drawstatus, 
         Plasma_storagedays, CSF_drawdate, CSF_drawage, CSF_drawstatus, CSF_storagedays, drawdate_diff, avg_drawdate) %>%
  as.data.frame()
  
# actually drawdate diff is all below 90!!! we can proceed without setting a time window. 
CSF_prot2_fil_z <- scale(CSF_protein2_fil[-c(1:3)])
plasma_prot2_fil_z <- scale(plasma_protein2_fil[-c(1:3)])
intake2_fil <- 1 - CSF_protein2_fil[-c(1:3)] / plasma_protein2_fil[-c(1:3)]

# we also need to change some diagnosis status 
# ADMCI -> AD, MCI-AD -> AD, HC -> CO

patientdata2_fil <- patientdata2_fil %>% 
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

patientdata2_AD_index <- patientdata2_fil$final_status == "AD" 
patientdata2_CO_index <- patientdata2_fil$final_status == "CO"
patientdata2_AD_index[is.na(patientdata2_AD_index)] <- FALSE
patientdata2_CO_index[is.na(patientdata2_CO_index)] <- FALSE
patientdata2_index <- patientdata2_AD_index | patientdata2_CO_index

patientdata2_fil <- patientdata2_fil[patientdata2_index, ]
intake2_fil <- intake2_fil[patientdata2_index, ]


# first, we proceed with differential analysis 
intake2_lm_models <- list()
for(i in 1:ncol(intake2_fil)) {
  intake2_lm_models[[i]] <- summary(lm(intake2_fil[, i] ~ 
                                         relevel(factor(patientdata2_fil$final_status), ref = "CO") + 
                                         relevel(factor(patientdata2_fil$Gender), ref = "F") +
                                         patientdata2_fil$avg_drawage + 
                                         patientdata2_fil$Plasma_storagedays + 
                                         patientdata2_fil$CSF_storagedays))
}

intake2_lm_summary <- data.frame()
for (i in 1:ncol(intake2_fil)) {
  intake2_tidy_result <- data.frame(intake2_lm_models[[i]]$coefficients) %>% 
    select(Estimate, Std..Error, Pr...t..) 
  colnames(intake2_tidy_result) <- c("Coefficient", "StdError", "Pval")
  intake2_tidy_result <- intake2_tidy_result[2:nrow(intake2_tidy_result), ] %>% 
    mutate(xVar = c("AD", "Male", "Age", "Plasma_storagedays", "CSF_storagedays"), 
           Protein = (colnames(intake2_fil))[i]) %>% 
    select(Protein, xVar, Coefficient, StdError, Pval)
  rownames(intake2_tidy_result) <- c()
  intake2_lm_summary <- rbind(intake2_lm_summary, intake2_tidy_result)
}


xVars2 <- c("AD", "Male", "Age", "Plasma_storagedays", "CSF_storagedays")

stanford_volcanodata <- list()
for (i in 1:length(xVars2)) {
  stanford_volcanodata[[i]] <- generate_volcanodata(intake2_lm_summary, variable = xVars2[i])
}
stanford_volcanoplot <- list()
for (i in 1:length(xVars2)) {
  stanford_volcanoplot[[i]] <- generate_volcanoplot(stanford_volcanodata[[i]], 
                                                          variable = xVars2[i])
}
gridExtra::grid.arrange(grobs = stanford_volcanoplot, ncol = 3)

for (i in 1:length(xVars2)) {
  write.csv(stanford_volcanodata[[i]], 
            paste0("data/generated_data/stanford_volcanodata/stanford_volcanodata_", xVars2[i]))
}











# visualization of intake 
generate_intakeplot <- function(protein, date_window, diagnosis) {
  if (missing(date_window)) {
    intake_lenient
  }
  if (missing(date_window)) {
    if (missing(diagnosis)) {
      
    } else if (diagnosis == "CO") {
      linedata <- data.frame(age = intake_age_seq, value = intake_pred_CO_z[, protein])
    } else if (diagnosis == "AD") {
      linedata <- data.frame(age = intake_age_seq, value = intake_pred_AD_z[, protein])
    }
  } else {
    if (date_window == "strict") {
      model <- loess()
    }
  }
  
  
  
  #choose correct loess interpolated 
  scatterdata <- intake_strict %>% select(protein)
  
  
  if (diagnosis == "AD" & lenient) {
    scatterdata <- intake_strict_AD %>% select(protein)
  } else if (diagnosis == "CO") {
    scatterdata <- intake_strict_CO %>% select(protein)
  }
  scatterdata <- intake_strict
  if (filter == "lenient" || filter == "strict") {
    scatterdata <- 
  }
  
  ggplot() + 
    geom_line(linedata, aes(x = age, y = value)) + 
    geom_point(scatterdata, aes(x = age, y = value))
}



# enricher: transcription factor focused analysis. 
# string-db: protein interaction network

# do same thing with stanford dataset too 
# 2 analysis: 1. clustering analysis over time / 
#             2. differential expression using linear models with sliding powers
# maybe accounting for different stages of disease in AD (measure of disease 
# severity independent of age)



knight_strict_volcanodata_outliers <- list()
for (i in 1:length(xVars)) {
  knight_strict_volcanodata_outliers[[i]] <- generate_volcanodata_outliers(intake_strict_lm_summary, 
                                                         variable = xVars[i])
}
knight_strict_volcanoplot_outliers <- list()
for (i in 1:length(xVars)) {
  knight_strict_volcanoplot_outliers[[i]] <- generate_volcanoplot(knight_strict_volcanodata_outliers[[i]], 
                                                         variable = xVars[i])
}

gridExtra::grid.arrange(grobs = knight_strict_volcanoplot_outliers, ncol = 3)



