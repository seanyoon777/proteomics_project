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
all_prots <- names(CSF_prots_fil[-c(1:3)])

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

# Create histogram of population 
drawdate_windows <- 15 * (0:6)
drawdate_hists <- list()
for (i in 1:6) {
  patientdata_fil_bin <- patientdata_fil[
    abs(patientdata_fil$drawdate_diff) <= drawdate_windows[i + 1] & 
      abs(patientdata_fil$drawdate_diff) >= drawdate_windows[i], ]
  
  drawdate_hists[[i]] <- patientdata_fil_bin %>% 
    ggplot(aes(x = avg_drawage, fill = final_status)) + 
    geom_histogram(binwidth = 1, breaks = seq(40, 100, by = 10), 
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
    scale_y_continuous(breaks = seq(0, 45, 5)) +
    scale_fill_manual(values = c("CO" = "#00BFC4", "AD" = "#F8766D"))
  
  maxcount <- max(hist(patientdata_fil_bin$avg_drawage, 
                       breaks = seq(30, 90, by = 10))$counts)
  drawdate_hists[[i]] <- drawdate_hists[[i]] + 
    annotate("text", x = 95, y = maxcount, label = paste0("n = ", nrow(patientdata_fil_bin)), 
             size = 3)
}
#AD = red, CO = blue
gridExtra::grid.arrange(grobs = drawdate_hists, ncol = 3)



# predict data and trends 
age_min <- max(min(patientdata_AD$avg_drawage), min(patientdata_CO$avg_drawage))
age_max <- min(max(patientdata_AD$avg_drawage), max(patientdata_CO$avg_drawage))
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

predratio_AD_z <- loess_predict(patientdata_AD$avg_drawage, ratio_AD_z, age_seq, "Age")
predratio_CO_z <- loess_predict(patientdata_CO$avg_drawage, ratio_CO_z, age_seq, "Age")

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
# we conclude both 10 clusters 

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
ratio_CO_z_clustnum <- proteinByClust(ratio_CO_z_clust, ratio_CO_z_dist, nclust)

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
    lmodels[[i]] <- summary(lm(protdata[, i] ~ patientdata$Gender + 
                                 patientdata$avg_drawage + 
                                 patientdata$Plasma_storagedays + 
                                 patientdata$CSF_storagedays))
  }
  return(lmodels)
}

ratio_AD_models <- generate_lmodels(ratio_AD, patientdata_AD)
ratio_CO_models <- generate_lmodels(ratio_CO, patientdata_CO)

generate_lmsummary <- function(lmodels) {
  lmsummary <- data.frame()
  for (i in 1:length(all_prots)) {
    tidyresult <- data.frame(lmodels[[i]]$coefficients) %>% 
      select(Estimate, Std..Error, Pr...t..)
    colnames(tidyresult) <- c("Coefficient", "StdError", "Pval")
    tidyresult <- tidyresult[2:nrow(tidyresult), ] %>% 
      mutate(xVar = c("Male", "Age", "Plasma_storagedays", "CSF_storagedays"), 
             Protein = all_prots[i]) %>% 
      dplyr::select(Protein, xVar, Coefficient, StdError, Pval)
    rownames(tidyresult) <- c()
    lmsummary <- rbind(lmsummary, tidyresult)
  }
  return(lmsummary)
}

ratio_AD_summary <- generate_lmsummary(ratio_AD_models)
ratio_CO_summary <- generate_lmsummary(ratio_CO_models)

xVars <- c("Male", "Age", "Plasma_storagedays", "CSF_storagedays")
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
generate_volcanoplot(ratio_AD_volcanodata_nonpadj)  # non-padj has pretty robust results

# Get enrichR pathways 



# knight before p-adjusted VS stanford volcano plots --> compare enrichR pathways 

#interaction model bw age and alzheimers disease status (extract CD? score)



