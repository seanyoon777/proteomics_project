#setting foundations
source("proteomics.R")
source("/Users/seonghyunyoon/Downloads/data_analysis_libraries/R_data_analysis.R")

rm(list = ls())

# 0. Basic data wrangling ----
## 0.1. Import data files ----
dir = "/Users/seonghyunyoon/Downloads/proteomics_project_1/"
barcodes <- read.csv(paste(dir, "CSF_plasma_matched_barcodes.csv", sep = "/"))

CSF <- read.csv(paste(dir, "CSF_metadata_2023-03-04.csv", sep = "/"))
CSF_prot <- read.csv(paste(dir, "CSFProts.log10.noLODFilter.csv", sep = "/"))
CSF_LOD <- read.csv(paste(dir, "CSF_ProteinMetadata_LOD_Samps_and_Dist_HumanHIVProtsOnly.csv", sep = "/"))

plasma <- read.csv(paste(dir, "Plasma_metadata_samplesOnly_2023-03-04.csv", sep = "/"))
plasma_prot <- read.csv(paste(dir, "plasmaProts.log10.csv", sep = "/"))
plasma_LOD <- read.csv(paste(dir, "ProteinMetadata_with_LOD.csv", sep = "/"))

## 0.3. Change CSF age (str) to numeric ----
CSF$Age <- gsub('^"&"$', '', CSF$Age)
CSF$Age <- as.numeric(CSF$Age)
  #raw_CSF$Age[is.na(CSF$Age)][1]
  #class(raw_CSF$Age[3])
  #raw_CSF$Age[3]

# 1. Plasma & CSF Relationship ----
## 1.1. Only leave genes present in both data ----
# plasma_red = only proteins that are present in CSF
# plasma_filtered = select rows with barcodes in barcodes$Plasma_Barcode
common_prots <- intersect(names(plasma_prot), names(CSF_prot))
CSF_red <- data.frame(Barcode = CSF$Barcode, CSF_prot[, common_prots])
plasma_red <- data.frame(Barcode = plasma_prot$Barcode, plasma_prot[, common_prots]) 

CSF_patientdata_filtered <- CSF[CSF$Barcode %in% barcodes$CSF_Barcode, ]
CSF_filtered <- CSF_red[CSF$Barcode %in% barcodes$CSF_Barcode, ]
plasma_filtered <- plasma_red %>% 
  filter(Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, barcodes$Plasma_Barcode)))

age_filtered <- plasma %>% 
  select(Age, Barcode) %>% 
  rename(Plasma_Barcode = Barcode) %>% 
  filter(Plasma_Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes$Plasma_Barcode)) %>% 
  mutate(CSF_Barcode = CSF_filtered$Barcode) %>% 
  select(Age, CSF_Barcode, Plasma_Barcode)

patientdata_filtered <- plasma %>% 
  rename(Plasma_Barcode = Barcode, Plasma_Zscore = ConnectivityZscore) %>% 
  filter(Plasma_Barcode %in% barcodes$Plasma_Barcode) %>% 
  arrange(match(Plasma_Barcode, barcodes$Plasma_Barcode)) %>% 
  mutate(CSF_Zscore = CSF_patientdata_filtered$ConnectivityZscore, 
         CSF_Barcode = CSF_patientdata_filtered$Barcode) %>% 
  .[, c(ncol(.), 1:(ncol(.) - 1))] %>% 
  .[, c(1:(ncol(.) - 2), ncol(.), ncol(.) - 1)]


## 1.2. Linear Regression on Intake to Predict Age ----
#PCA to CSF, PCA to plasma and then linear regression
PCA_CSF_filtered <- CSF_filtered %>% 
  select(-Barcode, -SampleId, -SampleType) %>% 
  prcomp(scale. = FALSE)

p_CSF_filtered_PCAimp <- fviz_eig(PCA_CSF_filtered) 
## 1.2.1. Clustering with Kmeans and with diagnosis group
PCA_CSF <- CSF_prot %>% 
  select(-SampleId, -SampleType) %>% 
  prcomp(scale. = FALSE)
PCA_CSF_2pc <- PCA_CSF$x[, 1:2]
PCA_CSF_clust <- kmeans(PCA_CSF_2pc, 4)
PCA_CSF_2pc <- data.frame(PCA_CSF_2pc, 
                                   group = as.factor(PCA_CSF_clust$cluster), 
                                   diagnosis_group = CSF$Diagnosis_group)

p_CSF_PCA <- PCA_CSF_2pc %>% 
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))

p_CSF_PCA_dg <- PCA_CSF_2pc %>% 
  ggplot(aes(x = PC1, y = PC2, color = diagnosis_group)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))

PCA_CSF_2pc %>% 
  subset(diagnosis_group %in% c("AD", "LBD", "MCI", "MCI-AD", "MCI-PD", "PD")) %>%
  ggplot(aes(x = PC1, y = PC2, color = diagnosis_group)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))

#TODO: heatmap with highest PC (col) & all the dep. variables that matter (row)----


PCA_CSF_filtered_2pc <- PCA_CSF_filtered$x[, 1:2]
plot(PCA_CSF_filtered_2pc)
PCA_CSF_filtered_clust <- kmeans(PCA_CSF_filtered_2pc, 4)
PCA_CSF_filtered_2pc <- data.frame(PCA_CSF_filtered_2pc, 
  group = as.factor(PCA_CSF_filtered_clust$cluster), 
  diagnosis_group = CSF_patientdata_filtered$Diagnosis_group)

p_CSF_reduced_PCA <- PCA_CSF_filtered_2pc %>% 
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))

p_CSF_reduced_PCA_dg <- PCA_CSF_filtered_2pc %>% 
  ggplot(aes(x = PC1, y = PC2, color = diagnosis_group)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))



PCA_plasma_filtered <- plasma_filtered %>% 
  select(-Barcode, -SampleId, -SampleType) %>% 
  prcomp(scale. = FALSE)

p_plasma_filtered_PCAimp <- fviz_eig(PCA_plasma_filtered) 

PCA_plasma_filtered_pc2 <- PCA_plasma_filtered$x[, 1:2]

LM_CSF_filtered <- lm(age_filtered$Age ~ PCA_CSF_filtered_pc2)
PCA_CSF_filtered_pc2 <- data.frame(Age = age_filtered$Age, PCA_CSF_filtered_pc2)
CSF_filtered_predage <- predict(LM_CSF_filtered, PCA_CSF_filtered_pc2[, 1:2])
p_CSF_filtered_predage_chronoage <- 
  ggplot(data = data.frame(Chronogical_Age = PCA_combined_filtered_pc[, 1], 
                           Predicted_Age = combined_filtered_predage), 
         aes(x = Chronogical_Age, y = Predicted_Age)) + 
  geom_point(color = "blue", size = 1, alpha = 0.2) + 
  geom_smooth(method = "lm", se = TRUE, color = "red", fullrange=TRUE) + 
  xlim(50, 80) +
  ylim(40, 100) + 
  ggtitle("Biological Age Prediction using CSF Proteomes (PCA)") + 
  xlab("Chronological Age (years)") + 
  ylab("Predicted Age (years)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))


# Try combined?
PCA_combined_filtered_pc <- data.frame(
  Age = age_filtered$Age, 
  PCA_CSF_filtered_pc2, 
  PCA_plasma_filtered_pc2)

LM_combined_filtered <- lm(Age ~ ., PCA_combined_filtered_pc)
combined_filtered_predage <- predict(LM_combined_filtered, PCA_combined_filtered_pc[, 2:ncol(PCA_combined_filtered_pc)])

p_combined_filtered_predage_chronoage <- 
  ggplot(data = data.frame(Chronogical_Age = PCA_combined_filtered_pc[, 1], 
                           Predicted_Age = combined_filtered_predage), 
                              aes(x = Chronogical_Age, y = Predicted_Age)) + 
  geom_point(color = "blue", size = 1, alpha = 0.2) + 
  geom_smooth(method = "lm", se = TRUE, color = "red", fullrange=TRUE) + 
  xlim(50, 80) +
  ylim(40, 100) + 
  ggtitle("Biological Age Prediction using Plasma and CSF Proteomes (PCA)") + 
  xlab("Chronological Age (years)") + 
  ylab("Predicted Age (years)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))


# No PCA since it doesn't really work, intake instead. 
intake_filtered <- (plasma_filtered[4:ncol(plasma_filtered)] - CSF_filtered[4:ncol(CSF_filtered)]) / plasma_filtered[4:ncol(plasma_filtered)]
PCA_intake_filtered <- prcomp(intake_filtered, scale. = FALSE)
p_intake_filtered_PCAimp <- fviz_eig(PCA_plasma_filtered) 

PCA_intake_filtered_2pc <- PCA_intake_filtered$x[, 1:2]

PCA_intake_filtered_2pc <- data.frame(
  Age = age_filtered$Age, 
  PC1 = PCA_intake_filtered_2pc[, 1], 
  PC2 = PCA_intake_filtered_2pc[, 2])

LM_intake_filtered <- lm(Age ~ ., PCA_intake_filtered_2pc)
intake_filtered_predage <- predict(LM_intake_filtered, PCA_intake_filtered_2pc[, 2:3])

p_intake_filtered_predage_chronoage <- 
  ggplot(data = data.frame(Chronogical_Age = PCA_intake_filtered_2pc[, 1], 
                           Predicted_Age = intake_filtered_predage), 
         aes(x = Chronogical_Age, y = Predicted_Age)) + 
  geom_point(color = "blue", size = 1, alpha = 0.2) + 
  geom_smooth(method = "lm", se = TRUE, color = "red", fullrange=TRUE) + 
  xlim(50, 80) +
  ylim(40, 100) + 
  ggtitle("Biological Age Prediction using Intake (PCA)") + 
  xlab("Chronological Age (years)") + 
  ylab("Predicted Age (years)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))


## 1.3. Use Intake for linear regression and create statistical testing table 
intake_lm_table <- vector("list", ncol(intake_filtered) - 3)
for(i in 1:(ncol(intake_filtered) - 3)) {
  intake_lm_table[[i]] <- summary(lm(intake_filtered[, i + 3] ~ 
                               patientdata_filtered$Diagnosis_group 
                             + patientdata_filtered$Age 
                             + patientdata_filtered$Gender
                             + patientdata_filtered$Storage_days 
                             + patientdata_filtered$Visit))
} # Very low r-squared values

#sigma

intake_lmer_table <- vector("list", ncol(intake_filtered) - 3)
for(i in 1:(ncol(intake_filtered) - 3)) {
  intake_lm_table[[i]] <- summary(lm(intake_filtered[, i + 3] ~ 
                                       patientdata_filtered$Diagnosis_group 
                                     + patientdata_filtered$Age 
                                     + patientdata_filtered$Gender
                                     + patientdata_filtered$Storage_days 
                                     + patientdata_filtered$Visit))
}

lmer(intake_filtered[, 3:ncol(intake_filtered)] ~ 
          patientdata_filtered$Diagnosis_group 
        + patientdata_filtered$Age 
        + patientdata_filtered$Gender
        + patientdata_filtered$Storage_days 
        + patientdata_filtered$Visit + (1 | PIDN))

## 1.4. Use Intake to create differential expressions ----
#lmer(prot_ratio ~ disease + age + gender + protein_storage + visit + (1|PIDN))



## 1.4. PCA and Clustering? (CSF, Plasma, Intake) ----
PCA_CSF_filtered_2pc <- PCA_CSF$x[, 1:2]
PCA_plasma_filtered_pc <- PCA_plasma_filtered$x[, 1:plasma_filtered_pca_ncomponents]

PCA_plasma_filtered <- 
  plasma_filtered %>%
  select(-Barcode, -SampleId, -SampleType) %>% 
  t() %>% prcomp()


## 1.4.1. Correlation Network (CSF) ----

CSF_cor <- CSF_prot %>% 
  select(-SampleId, -SampleType) %>% 
  correlate() %>%   # correlate() is equivalent to cor() but put NA as its diagonal entry and different class
  shave(upper = TRUE) %>%   # Shave the data frame to lower triangular matrix
  stretch(na.rm = TRUE) %>% 
  filter(r >= 0.99)

CSF_clust_elbow <- CSF_prot %>% 
  select(-SampleId, -SampleType) %>%
  fviz_nbclust(kmeans, method = "wss") +
  labs(subtitle = "Elbow method")

CSF_clust <- CSF_prot %>% 
  select(-SampleId, -SampleType) %>% t() %>% kmeans(4)

CSF_group <- CSF_prot %>% 
  select(-SampleId, -SampleType) %>% 
  names() %>% 
  data_frame(protein = ., group = as.factor(CSF_clust$cluster))

CSF_cor <- CSF_cor %>%
  as_tbl_graph(directed = FALSE) %>%
  activate(nodes) %>%
  left_join(CSF_group, by = c("name" = "protein")) %>%
  rename(label = name) %>% 
  activate(edges) %>% 
  rename(weight = r)

p_CSF_cornetwork <- ggraph(CSF_cor, layout = "fr") + 
  geom_edge_link(aes(width = weight), alpha = 0.2) + 
  scale_edge_width(range = c(0.2, 1.8)) +
  geom_node_point(aes(color = group), size = 2) + 
  geom_node_text(aes(label = label), size = 2, repel = TRUE, max.overlaps = 100) +
  theme_graph()

ggraph(CSF_cor, layout = "linear", circular = TRUE) +
  geom_node_point(aes(color = group), size = 2) +
  geom_edge_arc(aes(width = weight), alpha = 0.2) +
  scale_edge_width(range = c(0.2, 1.8)) +
  geom_node_text(aes(label = label), size = 2, repel = TRUE, max.overlaps = 100) +
  theme_graph()


## 1.4.2. Correlation Network (Plasma)
plasma_cor <- plasma_prot %>% 
  select(-Barcode, -SampleId, -SampleType) %>% 
  correlate() %>%   # correlate() is equivalent to cor() but                                                put NA as its diagonal entry and different class
  shave(upper = TRUE) %>%   # Shave the data frame to lower triangular matrix
  stretch(na.rm = TRUE) %>%           
  filter(r >= 0.99)

plasma_clust_elbow <- plasma_prot %>% 
  select(-Barcode, -SampleId, -SampleType) %>%
  fviz_nbclust(kmeans, method = "wss") +
  labs(subtitle = "Elbow method")

plasma_clust <- plasma_prot %>% 
  select(-Barcode, -SampleId, -SampleType) %>% t() %>% kmeans(3)

plasma_group <- plasma_prot %>% 
  select(-Barcode, -SampleId, -SampleType) %>% 
  names() %>% 
  data_frame(protein = ., group = as.factor(plasma_clust$cluster))

plasma_cor <- plasma_cor %>%
  as_tbl_graph(directed = FALSE) %>%
  activate(nodes) %>%
  left_join(plasma_group, by = c("name" = "protein")) %>%
  rename(label = name) %>% 
  activate(edges) %>% 
  rename(weight = r)

p_plasma_cornetwork <- ggraph(plasma_cor, layout = "fr") + 
  geom_edge_link(aes(width = weight), alpha = 0.2) + 
  scale_edge_width(range = c(0.2, 1.8)) +
  geom_node_point(aes(color = group), size = 2) + 
  geom_node_text(aes(label = label), size = 2, repel = TRUE, max.overlaps = 100) +
  theme_graph()

ggraph(plasma_cor, layout = "linear", circular = TRUE) +
  geom_node_point(aes(color = group), size = 2) +
  geom_edge_arc(aes(width = weight), alpha = 0.2) +
  scale_edge_width(range = c(0.2, 1.8)) +
  geom_node_text(aes(label = label), size = 2, repel = TRUE, max.overlaps = 100) +
  theme_graph()

# 2. Plasma & Age Relationship ----
## 2.1. Predicting age using plasma (PCA, lm) ----
  # Description: we use PCA to preprocess the features and then use linear 
  # regression. This is because the linear regression model on the raw dataset
  # is "not defined because of singularities" according to R, because the 
  # predictor variables are either linearly dependent or highly correlated with 
  # each other. Other models that are less sensitive to colinearity, such as LASSO 
  # or decision trees, have higher time complexity, which made them less suitable
  # to apply on a dataset like this, which has several thousand features. 
install.packages("pls")
library(pls)

# Process PCA model
plasma_chrono_age <- plasma$Age
plasma_prot <- data.frame(plasma_log10[4:ncol(plasma_log10)])
plasma_zscore <- scale(plasma_prot)
#plasma_prot_scaled <- plasma_prot %>%
#  preProcess(method = "scale") %>%
#  predict(plasma_prot)
plasma_pca_model <- prcomp(plasma_zscore, center = TRUE, scale. = TRUE)
plasma_pca_red <- data.frame(PC1 = plasma_pca_model$x[, 1], PC2 = plasma_pca_model$x[, 2])
p_plasma_pca <- ggplot(data = plasma_pca_red, aes(x = PC1, y = PC2)) +
  geom_point(size = 1, alpha = 0.2) + 
  ggtitle("PCA analysis of Plasma") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75), 
        plot.title = element_text(hjust = 0.5))

# Select num of clusters 
summary(plasma_pca_model) 
p_plasma_pca_importance <- fviz_eig(plasma_pca_model) 
#Q2. I got very high contribution from one variable ----
p_plasma_pca_clustnum <- fviz_nbclust(plasma_pca_model$x, FUNcluster=kmeans, k.max = 8) 
p_plasma_pca_clust <- eclust(plasma_pca_model$x, "kmeans", hc_metric = "eucliden", k = 2)

# Select num of components
plamsa_pca_var <- plasma_pca_model$sdev^2
plasma_pca_explained_var <- plamsa_pca_var / sum(plamsa_pca_var)
pca_cumulative_var <- cumsum(plasma_pca_explained_var)
plasma_age_ncomponents <- sum(pca_cumulative_var <= 0.95)



# Transform the predictor variables
plasma_pca_transformed <- predict(plasma_pca_model, newdata = plasma_zscore)[, 1:plasma_age_ncomponents]
plasma_age_lm <- train(Age ~ ., data = cbind(plasma_pca_transformed, Age = plasma$Age), method = "lm")
plasma_pred_age <- predict(plasma_age_lm, plasma_pca_transformed)

# Plot the predicted age and chronological age
summary(plasma_age_lm)
p_predage_chronoage <- ggplot(data = data.frame(Chronogical_Age = plasma_chrono_age, Predicted_Age = plasma_pred_age), 
       aes(x = Chronogical_Age, y = Predicted_Age)) + 
  geom_point(color = "blue", size = 1, alpha = 0.2) + 
  geom_smooth(method = "lm", se = TRUE, color = "red", fullrange=TRUE) + 
  xlim(40, 100) +
  ylim(40, 100) + 
  ggtitle("Prediction of Biological Age using Plasma Proteomes") + 
  xlab("Chronological Age (years)") + 
  ylab("Predicted Age (years)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))



#data exploration to see which features are highly correlated
plasma_prot_corr <- cor(plasma_prot)
n <- sum(plasma_prot_corr > 0.99)
idx <- which(plasma_prot_corr > 0.99, arr.ind = TRUE)

plasma_prot[1:10,1689]
ncol(plasma_prot)
cor_matclass(plasma_prot)
names(plasma_prot)[10:20]
tempidx <- grep("PNP", names(plasma_prot), value = TRUE)  #example of highly correlated 
plasma_prot[1:10, tempidx[1:2]]


#CHANGE OF PLASMA PROTEOMES BY AGE 
# merge the data frames
plasma_prot_zscore <- scale(plasma_prot)
plasma_prot_zscore_age <- data.frame(Age = plasma$Age, plasma_prot_zscore)

# convert protein columns to rows
plasma_long <- plasma_prot_zscore_age %>%
  gather(key = "protein", value = "zscore", -Age)

# plot the regression lines
p_plasma_age_zscore <- plasma_long %>% 
  ggplot(aes(x = Age, y = zscore, group = protein)) + 
  geom_line(stat = "smooth", method = "loess", se = FALSE, color = "black", 
            size = 0.25, alpha = 0.3) + 
  ggtitle("Protein profiles (Plasma)") +
  xlab("Age (years)") + 
  ylab("Protein Level (z-score)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))
  
#CHANGE OF CSF PROTEOMES BY AGE 
# merge the data frames
CSF_prot <- CSF_log10_noLOD[, 3:ncol(CSF_log10_noLOD)]
CSF_prot_zscore <- scale(CSF_prot)
CSF_prot_zscore_age <- data.frame(Age = CSF$Age, CSF_prot_zscore)

# convert protein columns to rows
CSF_long <- CSF_prot_zscore_age %>%
  gather(key = "protein", value = "zscore", -Age)

# plot the regression lines
p_CSF_age_zscore <- CSF_long %>% 
  ggplot(aes(x = Age, y = zscore, group = protein)) + 
  geom_line(stat = "smooth", method = "loess", se = FALSE, color = "black", 
            size = 0.25, alpha = 0.4) + 
  ggtitle("Protein profiles (CSF)") +
  xlab("Age (years)") + 
  ylab("Protein Level (z-score)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        plot.title = element_text(hjust = 0.5))






