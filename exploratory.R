#setting foundations
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
install.packages("corrplot")

library(tidyverse)
library(purrr)
library(ggplot2)
library(corrplot)
library(caret)
library(dplyr)
library(data.table)
rm(list = ls())

#import data files
dir = "/Users/seonghyunyoon/Downloads/proteomics_project_1/"
barcodes <- read.csv(paste(dir, "CSF_plasma_matched_barcodes.csv", sep = "/"))

raw_CSF <- read.csv(paste(dir, "CSF_metadata_2023-03-04.csv", sep = "/"))
raw_CSF_log10_noLOD <- read.csv(paste(dir, "CSFProts.log10.noLODFilter.csv", sep = "/"))
raw_CSF_LOD <- read.csv(paste(dir, "CSF_ProteinMetadata_LOD_Samps_and_Dist_HumanHIVProtsOnly.csv", sep = "/"))

raw_plasma <- read.csv(paste(dir, "Plasma_metadata_samplesOnly_2023-03-04.csv", sep = "/"))
raw_plasma_log10 <- read.csv(paste(dir, "plasmaProts.log10.csv", sep = "/"))
raw_plasma_LOD <- read.csv(paste(dir, "ProteinMetadata_with_LOD.csv", sep = "/"))

#change order
CSF_order <- order(match(raw_CSF$Barcode, barcodes$CSF_Barcode))
CSF <- raw_CSF[CSF_order, ]
CSF_log10_noLOD <- raw_CSF_log10_noLOD[CSF_order, ]

plasma_order <- order(match(raw_plasma$Barcode, barcodes$Plasma_Barcode))
plasma <- raw_plasma[plasma_order, ]
plasma_log10 <- raw_plasma_log10[plasma_order, ]

CSF$Age <- gsub('^"&"$', '', CSF$Age)
CSF$Age <- as.numeric(CSF$Age)
raw_CSF$Age[is.na(CSF$Age)][1]
class(raw_CSF$Age[3])
raw_CSF$Age[3]

#ANALYZING RELATIONSHIP B/W PLASMA AND CSF COMPOSITION
#remove genes not present in both data
common_cols <- intersect(names(plasma_log10), names(CSF_log10_noLOD))
plasma_red <- plasma_log10[, common_cols]
CSF_red <- CSF_log10_noLOD[, common_cols]

# PREDICTING AGE USING PLASMA PROTEOMES (PCA PREPROCESSING)
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
plasma_prot <- plasma_log10[4:ncol(plasma_log10)]
plasma_prot_scaled <- plasma_prot %>%
  preProcess(method = "scale") %>%
  predict(plasma_prot)
plasma_pca_model <- prcomp(plasma_prot_scaled, center = TRUE, scale. = TRUE)

# Select num of components
plamsa_pca_var <- plasma_pca_model$sdev^2
plasma_pca_explained_var <- plamsa_pca_var / sum(plamsa_pca_var)
pca_cumulative_var <- cumsum(plasma_pca_explained_var)
plasma_age_ncomponents <- sum(pca_cumulative_var <= 0.95)

# Transform the predictor variables
plasma_pca_transformed <- predict(plasma_pca_model, newdata = plasma_prot_scaled)[, 1:plasma_age_ncomponents]
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
tempidx <- grep("PNP", names(plasma_prot), value = TRUE)#example of highly correlated 
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








data(iris)
x <- iris[, 1:4]


# Calculate correlation matrix
corr <- cor(x)

# Convert to adjacency matrix
adj <- as.matrix(corr)

# Create network graph
install.packages("igraph")
library(igraph)
g <- graph.adjacency(adj, mode = "undirected", weighted = TRUE)
hist(x = CSF$Age)
head(plasma$Age)
# Plot network graph
plot(g, layout = layout.fruchterman.reingold, vertex.label = colnames(x))


