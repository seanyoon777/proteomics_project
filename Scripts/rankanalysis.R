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
all_patientdata <- get_biodata("combined/patientdata.csv")
all_patientdata_AD <- get_biodata("combined/patientdata_AD.csv")
all_patientdata_CO <- get_biodata("combined/patientdata_CO.csv")

all_plasma <- get_biodata("combined/plasma_raw_lodfiltered.csv")
all_plasma_AD <- get_biodata("combined/plasma_raw_AD_lodfiltered.csv")
all_plasma_CO <- get_biodata("combined/plasma_raw_CO_lodfiltered.csv")

all_CSF <- get_biodata("combined/CSF_raw_lodfiltered.csv")
all_CSF_AD <- get_biodata("combined/CSF_raw_AD_lodfiltered.csv")
all_CSF_CO <- get_biodata("combined/CSF_raw_CO_lodfiltered.csv")

# Experiment 1: Raw values  
all_CSF_plasma <- all_CSF / all_plasma
all_CSF_plasma_AD <- all_CSF_AD / all_plasma_AD
all_CSF_plasma_CO <- all_CSF_CO / all_plasma_CO

# Experiment 2 default: z-score, rank norm
rank_transform <- function(df) {
  temp <- as.data.frame(scale(df))
  temp <- as.data.frame(t(apply(temp, 1, rank)))
  return(temp)
}

all_plasma <- rank_transform(all_plasma)
all_plasma_AD <- rank_transform(all_plasma_AD)
all_plasma_CO <- rank_transform(all_plasma_CO)
all_CSF <- rank_transform(all_CSF)
all_CSF_AD <- rank_transform(all_CSF_AD)
all_CSF_CO <- rank_transform(all_CSF_CO)

# Experiment 2a: rank ratio
all_CSF_plasma <- all_CSF / all_plasma
all_CSF_plasma_AD <- all_CSF_AD / all_plasma_AD
all_CSF_plasma_CO <- all_CSF_CO / all_plasma_CO

# Experiment 2b: rank difference
all_CSF_plasma <- all_CSF - all_plasma
all_CSF_plasma_AD <- all_CSF_AD - all_plasma_AD
all_CSF_plasma_CO <- all_CSF_CO - all_plasma_CO



rowwise_inverse_normal <- function(row) {
  n <- length(row)
  ranks <- rank(row, ties.method = "average")
  percentiles <- (ranks - 0.5) / n
  return(qnorm(percentiles))
}


# all_plasma_rank <- as.data.frame(t(apply(all_plasma_scale, 1, rowwise_inverse_normal)))
# all_CSF_rank <- as.data.frame(t(apply(all_CSF_scale, 1, rowwise_inverse_normal)))



# Define DEA function
dea <- function(protdata, patientdata, xVars, xLabs_lmsummary, xLabs_plot, 
                ncol_plot, prot_num, reverse_bools) {
  lmodels <- generate_lmodels(protdata, patientdata, xVars)
  lsummary <- generate_lmsummary(lmodels, colnames(protdata), xLabs_lmsummary)
  volcanodata <- generate_volcanodata(lsummary, xLabs_plot)
  generate_volcanoplot(volcanodata, xLabs_plot, ncol_plot, prot_num, reverse_bools)
}


# AD/CO included 
xVars <- c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect")
xLabs_lmsummary <- c("Male", "Age", "AD", "Drawdate difference", "storage days", "Study bias")
xLabs_plot <- c("Male", "Age", "AD")
reverse_bools <- c(FALSE, FALSE, TRUE)
dedata_ratio_all <- dea(all_CSF_plasma, all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 16, reverse_bools)

# AD/CO Separated Volcanoplots
xVars <- c("Sex", "avg_drawage", "drawdate_diff", "storage_days", "batch_effect")
xLabs_lmsummary <- c("Male", "Age", "Drawdate difference", "Storage days", "Study bias")
xLabs_plot <- c("Male", "Age")
reverse_bools <- c(FALSE, FALSE)
dedata_ratio_AD <- dea(all_CSF_plasma_AD, all_patientdata_AD, xVars, xLabs_lmsummary, xLabs_plot, 2, 16, reverse_bools)
dedata_ratio_CO <- dea(all_CSF_plasma_CO, all_patientdata_CO, xVars, xLabs_lmsummary, xLabs_plot, 2, 16, reverse_bools)

# CSF & Plasma Volcano plots 
xVars <- c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect")
xLabs_lmsummary <- c("Male", "Age", "AD", "drawdate_diff", "storage_days", "batch_effect")
xLabs_plot <- c("Male", "Age", "AD")
reverse_bools <- c(FALSE, FALSE, TRUE)
dedata_CSF_all <- dea(all_CSF, all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 16, reverse_bools)
dedata_plasma_all <- dea(all_plasma, all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 16, reverse_bools)

# Interaction model 


# Robustness Analysis: Volcanoplots for 20 & 40 day interval drawdate bins
drawdate_bins_20 <- c(0, 1, 20, 40, 60, 80, 100, 120)
dedata_by_bin(all_CSF_plasma, all_patientdata, drawdate_bins_40, 1, ncol = 2)

protdata, patientdata, drawdate_bins, factor_num, ncol, prot_num = 10)

dea_by_bin(all_CSF, all_patientdata, drawdate_bins_20, 1, ncol = 4)
drawdate_bins_40 <- c(0, 1, 40, 80, 120)
dea_by_bin(all_CSF_plasma, all_patientdata, drawdate_bins_40, 1, ncol = 2)

# Robustness Analysis: log2fc vs log2fc for 20 & 40 day interval drawdate bins


# Robustness Analysis: p value vs p value for 20 & 40 day interval drawdate bins





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


