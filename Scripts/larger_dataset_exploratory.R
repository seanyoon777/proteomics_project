# 0. Load data ----
source("proteomics.R")

rm(list = ls())

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

CSF_meta <- match_by(CSF_meta, CSF_patient)
CSF_prot <- match_by(CSF_prot, CSF_patient)
CSF_prot_zscored <- CSF_prot
CSF_prot_zscored[all_prots] <- scale(CSF_prot[all_prots])

## 0.2. Wrangle intake dataset ----
common_ID <- intersect(CSF_patient$PA_DB_UID, plasma_patient$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .)

CSF_meta_fil <- match_by(CSF_meta, common_ID)
CSF_prot_fil <- match_by(CSF_prot, common_ID)
CSF_patient_fil <- match_by(CSF_patient, common_ID)
plasma_meta_fil <- match_by(plasma_meta, common_ID)
plasma_prot_fil <- match_by(plasma_prot, common_ID)
plasma_patient_fil <- match_by(plasma_patient, common_ID)
patient_meta_fil <- match_by(patient_meta, common_ID)

intake_fil <- 1 - CSF_prot_fil[4:nprots] / plasma_prot_fil[4:nprots]

## 0.2. Wrangle common patientdata
drawdate_diff <- difftime(CSF_patient_fil$DrawDate, plasma_patient_fil$drawdate, units = "days") %>%
  data.frame(Diff = ., absDiff = abs(.))

# draw histograms of drawdate window
hist(as.numeric(drawdate_diff$absDiff))
drawdate_max <- max(as.numeric(drawdate_diff$absDiff))
hist(as.numeric(drawdate_diff$absDiff), breaks = 30 * (0:(drawdate_max / 30 + 1)))  # smaller bins
drawdate_smallerthan180 <- as.numeric(drawdate_diff$absDiff) %>% # zoomed in, until 6 months
  data.frame(diff = .) %>% filter(diff <= 180) 
  hist(as.numeric(drawdate_smallerthan180$diff), breaks = 15 * (0:12)) 

# summary of data
drawdate_windows <- c(0, 30 * (1:6))
summary_drawdate_diff <- sapply(drawdate_windows, 
                                function(x) sum(as.numeric(drawdate_diff$absDiff) <= x)) %>% 
  data.frame(drawdate_windows, frequency = .)

patientdata_fil <- merge(CSF_patient_fil, plasma_patient_fil, by = "PA_DB_UID") %>% 
  filter(Last.status %in% c("AD", "CO")) %>% 
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status) %>%
  mutate(drawdate_diff = abs(difftime(Plasma_drawdate, CSF_drawdate, units = "days")))


# TO DO: visualizations of distribution of ages at different cutoffs of drawdate windows
# TO DO: AD / CO too (ridge or stacked histogram?) 
patientdata_fil_summary <- patientdata_fil[patientdata_fil$Plasma_drawstatus == 
                                             patientdata_fil$CSF_drawstatus, ] %>% 
  mutate(drawage = (Plasma_drawage + CSF_drawage) / 2)
  

# Create a list to store the plots
drawdate_windows <- 30 * (0:12)
patientdata_fil_hists <- list()
for (i in 1:12) {
  patientdata_fil_summary_bin <- patientdata_fil_summary[
    patientdata_fil_summary$drawdate_diff <= drawdate_windows[i + 1] & 
      patientdata_fil_summary$drawdate_diff >= drawdate_windows[i], ]
  patientdata_fil_hists[[i]] <- patientdata_fil_summary_bin %>% 
   ggplot(aes(x = drawage, fill = Plasma_drawstatus)) + 
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
    scale_y_continuous(breaks = seq(0, 45, 5)) 
  maxcount <- max(hist(patientdata_fil_summary_bin$drawage, 
                       breaks = seq(40, 100, by = 10))$counts)
  patientdata_fil_hists[[i]] <- patientdata_fil_hists[[i]] + 
    annotate("text", x = 95, y = maxcount, label = paste0("n = ", nrow(patientdata_fil_summary_bin)), 
             size = 3)
}
#AD = red, CO = blue
gridExtra::grid.arrange(grobs = patientdata_fil_hists, ncol = 4)


# 1. LME model----
# 1.1. LME of intake & Visualization


# Zero intercept model --> instead of picking specific reference factor, 
# sets mean value of all data as intercept. --> look into this 
# Run LM and / LMER if there's enough different number of visits 
# Visualize tables: volcano plot with coffecients in x axis and -log10(q) 
# in y axis (preferably w/ labels and w/ y threshold value) 


# 2. CSF AD/CO Clustering & Heatmap
# Separating CO / AD data
split_CSF_patient <- split(CSF_patient, CSF_patient$Last.status)
CSF_age_CO <- split_CSF_patient$CO[, 5]
CSF_age_AD <- split_CSF_patient$AD[, 5]
split_CSF_prot_z <- split(CSF_prot_zscored, CSF_patient$Last.status)
CSF_prot_z_CO <- split_CSF_prot_z$CO[, 4:ncol(split_CSF_prot_z$CO)]
CSF_prot_z_AD <- split_CSF_prot_z$AD[, 4:ncol(split_CSF_prot_z$AD)]
CSF_max_age <- min(max(CSF_age_CO), max(CSF_age_AD))
CSF_min_age <- max(min(CSF_age_CO), min(CSF_age_AD))
CSF_age_resolution <- 0.25
CSF_age_seq <- seq(CSF_min_age, CSF_max_age, by = CSF_age_resolution)
CSF_AD_age_seq <- seq(min(CSF_age_AD), max(CSF_age_AD), by = CSF_age_resolution)

# Clustering proteins for AD 
CSF_prot_AD_models <- vector("list", length(all_prots))
for (i in 1:length(all_prots)) {
  CSF_prot_AD_models[[i]] <- loess(CSF_prot_z_AD[, i] ~ CSF_age_AD)
}
CSF_prot_AD_preddata <- data.frame(axis.values.x = CSF_AD_age_seq)
for (i in 1:length(all_prots)) {
  CSF_prot_AD_preddata[, i + 1] <- predict(CSF_prot_AD_models[[i]], CSF_AD_age_seq)
}
names(CSF_prot_AD_preddata) <- c("axis.values.x", all_prots)

CSF_prot_CO_models <- vector("list", length(all_prots))
for (i in 1:length(all_prots)) {
  CSF_prot_CO_models[[i]] <- loess(CSF_prot_z_CO[, i] ~ CSF_age_CO)
}
CSF_prot_CO_preddata <- data.frame(axis.values.x = CSF_age_seq)
for (i in 1:length(all_prots)) {
  CSF_prot_CO_preddata[, i + 1] <- predict(CSF_prot_CO_models[[i]], CSF_age_seq)
}
names(CSF_prot_CO_preddata) <- c("axis.values.x", all_prots)


CSF_silhouette_list <- numeric(14)
CSF_prot_AD_preddata_dist <- dist(CSF_prot_AD_preddata)
CSF_prot_AD_preddata_clust <- hclust(CSF_prot_AD_preddata_dist)
for (k in 2:15) {
  CSF_cut_clusters <- cutree(CSF_prot_AD_preddata_clust, k)
  CSF_vals <- silhouette(CSF_cut_clusters, CSF_prot_AD_preddata_dist)
  CSF_silhouette_list[k - 1] <- mean(CSF_vals[, 3])
}
plot(CSF_silhouette_list) # optimal number of clusters = 8
CSF_nclust <- 8

CSF_prot_AD_clust <- cutree(CSF_prot_AD_preddata_clust, k = CSF_nclust) %>%
  data.frame(cluster = .)
CSF_prot_AD_clust$variables <- all_prots
CSF_prot_AD_long <- CSF_prot_AD_preddata %>% gather(key = variables, value = value, -axis.values.x)
CSF_prot_AD_clust <- CSF_prot_AD_long %>% inner_join(CSF_prot_AD_clust, by = "variables") %>%
  arrange(cluster, variables, axis.values.x)

# Create heatmap
CSF_plot_agebin <- as.character(CSF_age_seq)
CSF_plot_agebin[!round(as.numeric(CSF_plot_agebin), 1) %% 10 == 0] <- ""

CSF_protreordered_AD_preddata <- CSF_prot_AD_preddata


CSF_prot_col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("lightblue", "black", "Yellow"))
CSF_prot_AD_heatdata <- t(as.matrix(CSF_prot_AD_preddata[-1]))
colnames(CSF_prot_AD_heatdata) <- CSF_plot_agebin
CSF_AD_heatmap <- Heatmap(CSF_prot_AD_heatdata, col = CSF_prot_col_fun, 
                          show_row_names = FALSE, show_column_dend = FALSE,
                          row_dend_reorder = TRUE, 
                          show_row_dend = FALSE, 
                          column_order = 1:ncol(CSF_prot_AD_heatdata), name = "Z-score", 
                          column_title = "CSF Protein Trajectories (AD)", 
                          row_title = "7,596 CSF Proteins", 
                          heatmap_legend_param = list(direction = "horizontal"))
CSF_prot_AD_order <- row_order(CSF_AD_heatmap)
CSF_prot_AD_heatdata <- CSF_prot_AD_heatdata[, CSF_AD_age_seq >= CSF_min_age & 
                                         CSF_AD_age_seq <= CSF_max_age] 
CSF_AD_heatmap <- 
  Heatmap(CSF_prot_AD_heatdata, col = CSF_prot_col_fun, show_row_names = FALSE, 
          show_column_dend = FALSE,
          row_order = CSF_prot_AD_order, show_row_dend = FALSE,
          column_order = 1:ncol(CSF_prot_AD_heatdata), show_heatmap_legend = FALSE, 
          column_title = "Age (Years)", row_title = "7,596 CSF Proteins", 
          column_title_side = "bottom", 
          column_names_rot = 0)
CSF_heatmap_legend <- Legend(col_fun = CSF_prot_col_fun, title = "Z-score", 
                             at = c(-0.5, 0, 0.5),
                             direction = "horizontal", title_position = "topcenter")
draw(CSF_AD_heatmap, 
     column_title= "CSF Protein Trajectories (AD)",
     column_title_gp=grid::gpar(fontsize=16), 
     heatmap_legend_list = list(CSF_heatmap_legend), heatmap_legend_side = "bottom")

CSF_prot_CO_heatdata <- t(as.matrix(CSF_prot_CO_preddata[-1]))
colnames(CSF_prot_CO_heatdata) <- CSF_plot_agebin
CSF_CO_heatmap <- Heatmap(CSF_prot_CO_heatdata, col = CSF_prot_col_fun, show_row_names = FALSE, 
                          show_column_dend = FALSE,
                          row_order = CSF_prot_AD_order, show_row_dend = FALSE,
                          column_order = 1:ncol(CSF_prot_AD_heatdata), show_heatmap_legend = FALSE, 
                          column_title = "Age (Years)", row_title = "7,596 CSF Proteins", 
                          column_title_side = "bottom", 
                          column_names_rot = 0)
draw(CSF_CO_heatmap, 
     column_title= "CSF Protein Trajectories (CO)",
     column_title_gp=grid::gpar(fontsize=16), 
     heatmap_legend_list = list(CSF_heatmap_legend), heatmap_legend_side = "bottom")

CSF_CO_heatmap <- Heatmap(CSF_prot_CO_heatdata, col = CSF_prot_col_fun, show_row_names = FALSE, 
                          show_column_dend = FALSE,
                          row_dend_reorder = TRUE, show_row_dend = FALSE,
                          column_order = 1:ncol(CSF_prot_AD_heatdata), show_heatmap_legend = FALSE, 
                          column_title = "Age (Years)", row_title = "7,596 CSF Proteins", 
                          column_title_side = "bottom", 
                          column_names_rot = 0)
draw(CSF_CO_heatmap, 
     column_title= "CSF Protein Trajectories (CO)",
     column_title_gp=grid::gpar(fontsize=16), 
     heatmap_legend_list = list(CSF_heatmap_legend), heatmap_legend_side = "bottom")



# Test clustering
Figure_CSF_prot_AD_clust <- CSF_prot_AD_clust %>%
  ggplot(aes(x = axis.values.x, y = value, group = variables, color = factor(cluster))) +
  geom_line(stat = "smooth", method = "loess", se = FALSE, linewidth = 0.5, alpha = 0.1) +
  scale_color_manual(
    values = c("#FF6F61", "#6B5B95", "#88B04B", "#FFA500", "#92A8D1", "#FF69B4", "#955251")) +
  labs(x = "X-axis", y = "Y-axis", color = "Cluster") +
  facet_wrap(~ cluster, ncol = 3, nrow = 4, labeller = labeller(cluster = as.character)) +
  stat_summary(fun.data = "mean_cl_normal", geom = "line",
               aes(group = cluster, color = factor(cluster)),
               linewidth = 0.8, alpha = 1, linetype = "solid", color = "white") +
  stat_summary(fun.data = "mean_cl_normal", geom = "line",
               aes(group = cluster, color = factor(cluster)),
               linewidth = 0.5, alpha = 0.8)












