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
drawdate_windows <- c(1, 30 * (1:6))
summary_drawdate_diff <- sapply(drawdate_windows, 
                                function(x) sum(as.numeric(drawdate_diff$absDiff) <= x)) %>% 
  data.frame(drawdate_windows, frequency = .)

patientdata_fil <- merge(CSF_patient_fil, plasma_patient_fil, by = "PA_DB_UID") %>% 
  filter(Last.status %in% c("AD", "CO")) %>% 
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status)


# TO DO: visualizations of distribution of ages at different cutoffs of drawdate windows
# TO DO: AD / CO too (ridge or stacked histogram?)

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

# Clustering proteins for AD, Create heatmap
CSF_prot_AD_preddata <- data.loess(CSF_age_AD, CSF_prot_z_AD)
CSF_prot_AD_preddata_dist <- dist(t(CSF_prot_AD_preddata[-1]), method="euclidean")
CSF_prot_AD_preddata_clust <- hclust(CSF_prot_AD_preddata_dist, method = "ward.D")

CSF_silhouette_list <- numeric(9)

for (k in 2:10) {
  CSF_cut_clusters <- cutree(CSF_prot_AD_preddata_clust, k)
  CSF_vals <- silhouette(CSF_cut_clusters, CSF_prot_AD_preddata_dist)
  CSF_silhouette_list[k - 1] <- mean(CSF_vals[, 3])
}
plot(CSF_silhouette_list) # optimal number of clusters = 5
CSF_nclust <- 5

CSF_prot_AD_clust <- cutree(CSF_prot_AD_preddata_clust, k = CSF_nclust) %>%
  data.frame(cluster = .)
CSF_prot_AD_clust$variables <- all_prots
CSF_prot_AD_long <- CSF_prot_AD_preddata %>% gather(key = variables, value = value, -axis.values.x)
CSF_prot_AD_clust <- CSF_prot_AD_long %>% inner_join(CSF_prot_AD_clust, by = "variables") %>%
  arrange(cluster, variables, axis.values.x)


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

# Create heatmap
CSF_prot_AD_order <- CSF_prot_AD_preddata_clust$order
CSF_protreordered_AD_preddata <- CSF_prot_AD_preddata
CSF_protreordered_AD_preddata[-1] <- select(CSF_prot_AD_preddata[-1], CSF_prot_AD_order) 
CSF_protreordered_AD_preddata_long <- CSF_protreordered_AD_preddata %>% 
  gather(key = variables, value = value, -axis.values.x)


heatmap_colorpalette <- CSF_prot_col_fun(seq(0.5, -0.5, by=-0.01))



CSF_prot_col_fun <- function(x, range) {
  if (x <= range[1]) {
    return(range[1])
  } else if (x >= range[2]) {
    return(range[2])
  }
  return(x)
}

AD_CSF_heatmap <- CSF_protreordered_AD_preddata_long %>% 
  mutate(value = CSF_prot_col_fun(value, c(-1, 1))) %>% 
  ggplot(aes(x = axis.values.x, y = variables, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank())

AD_CSF_heatmap

scale_fill_gradient(colors = my_color(100), na.value = "grey90") 


ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)

CSF_prot_col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
CSF_prot_AD_heatdata <- CSF_prot_AD_preddata[-1]
CSF_prot_AD_heatdata <- t(as.matrix(CSF_prot_AD_heatdata[CSF_prot_AD_order, ]))
  
CSF_AD_heatmap <- Heatmap(CSF_prot_AD_heatdata, col = CSF_prot_col_fun, 
        show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = TRUE)

CSF_AD_heatmap

ggsave("CSF_AD_heatmap.png", CSF_AD_heatmap, width = 17, height = 9.6)

Heatmap(mat, name = "mat", row_order = sort(rownames(mat)), 
        column_order = sort(colnames(mat)),
        column_title = "reorder matrix by row/column names")





ht_list <- Heatmap(as.matrix(CSF_prot_AD_preddata[-1]), name = "Age", col = CSF_prot_col_fun, 
                   column_title = "CSF Protein Expression Level") 


draw(ht_list)



  Heatmap(direction, name = "direction", col = direction_col) +
  Heatmap(mat_expr[, column_tree$order], name = "expression", 
          col = expr_col_fun, 
          column_order = column_order, 
          top_annotation = ha2, column_title = "Expression") +
  Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
  Heatmap(gene_type, name = "gene type", col = gene_type_col) +
  Heatmap(anno_gene, name = "anno_gene", col = anno_gene_col) +
  Heatmap(dist, name = "dist_tss", col = dist_col_fun) +
  Heatmap(anno_enhancer, name = "anno_enhancer", col = enhancer_col_fun, 
          cluster_columns = FALSE, column_title = "Enhancer")

draw(ht_list, row_km = 2, row_split = direction,
     column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "bottom")


# Same clusters for CO, Create heatmap


