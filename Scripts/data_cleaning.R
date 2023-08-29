
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
  mutate(CSF_finaldate = BirthYR + age.at.last, plasma_finaldate = birth_Year + Age_at_last) %>% 
  mutate(CSF_finaldiff = CSF_finaldate - as.numeric(format(DrawDate, format="%Y")), 
         plasma_finaldiff = plasma_finaldate - as.numeric(format(drawdate, format="%Y")))  %>%
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         CSF_drawstatus = CC_at_LP_draw, CSF_finalstatus = Last.status, 
         plasma_drawstatus = CDR_status_at_blood_draw, plasma_finalstatus = final_cc_status.updated,
         CSF_drawdate = DrawDate, CSF_finaldate, 
         plasma_drawdate = drawdate, plasma_finaldate, CSF_finaldiff, plasma_finaldiff,
         CSF_drawage = age_at_csf_draw, plasma_drawage = Age_at_blood_draw.updated
         ) %>%
  mutate(drawdate_diff = abs(difftime(plasma_drawdate, CSF_drawdate, units = "days")), 
         avg_drawdate = as.Date(mean(c(plasma_drawdate, CSF_drawdate))))


## PCA Analysis
knight_index_120 <- knight_patients$drawdate_diff <= 120
knight_patients <- knight_patients[knight_index_120, ]
knight_CSF_prots <- knight_CSF_prots[knight_index_120, ]
knight_plasma_prots <- knight_plasma_prots[knight_index_120, ]

all_prots <- intersect(names(knight_plasma_prots)[4:ncol(knight_plasma_prots)], 
                       names(knight_CSF_prots)[4:ncol(knight_CSF_prots)])
nprots <- length(all_prots)
  
knight_ratios <- knight_CSF_prots[-c(1:3)] / knight_plasma_prots[-c(1:3)] 

knight_confirmed_CO <- knight_patients$CSF_drawstatus == "CO" & knight_patients$CSF_finalstatus == "CO" &
  knight_patients$plasma_drawstatus == "CO" & knight_patients$plasma_finalstatus == "CO" 
knight_confirmed_AD <- knight_patients$CSF_drawstatus == "CA" & knight_patients$CSF_finalstatus == "AD" &
  knight_patients$plasma_drawstatus == "AD" & knight_patients$plasma_finalstatus == "AD" 
knight_confirmed <- knight_confirmed_CO | knight_confirmed_AD

knight_patients_confirmed <- knight_patients[knight_confirmed, ]
knight_plasma_confirmed <- knight_plasma_prots[knight_confirmed, ]
knight_CSF_confirmed <- knight_CSF_prots[knight_confirmed, ]
knight_ratios_confirmed <- knight_ratios[knight_confirmed, ]

knight_ratios_confirmed_pca <- prcomp(knight_ratios_confirmed, center = TRUE, scale. = TRUE)
summary(knight_ratios_confirmed_pca)

knight_ratios_confirmed_pca_scores <- as.data.frame(knight_ratios_confirmed_pca$x)
knight_ratios_confirmed_pca_combined <- cbind(knight_ratios_confirmed_pca_scores, 
                                              Species = knight_patients_confirmed$CSF_finalstatus)
ggplot(knight_ratios_confirmed_pca_combined, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Ratio",
       x = "Principal Component 1",
       y = "Principal Component 2")

knight_plasma_confirmed_pca <- prcomp(knight_plasma_confirmed[-c(1:3)], center = TRUE, scale. = TRUE)
summary(knight_plasma_confirmed_pca)
knight_plasma_confirmed_pca_scores <- as.data.frame(knight_plasma_confirmed_pca$x)
knight_plasma_confirmed_pca_combined <- cbind(knight_plasma_confirmed_pca_scores, 
                                              Species = knight_patients_confirmed$CSF_finalstatus)
ggplot(knight_plasma_confirmed_pca_combined, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Plasma",
       x = "Principal Component 1",
       y = "Principal Component 2")

knight_CSF_confirmed_pca <- prcomp(knight_CSF_confirmed[-c(1:3)], center = TRUE, scale. = TRUE)
summary(knight_CSF_confirmed_pca)
knight_CSF_confirmed_pca_scores <- as.data.frame(knight_CSF_confirmed_pca$x)
knight_CSF_confirmed_pca_combined <- cbind(knight_CSF_confirmed_pca_scores, 
                                              Species = knight_patients_confirmed$CSF_finalstatus)
ggplot(knight_CSF_confirmed_pca_combined, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of CSF",
       x = "Principal Component 1",
       y = "Principal Component 2")

knight_CSF_confirmed_pca_var <- (knight_CSF_confirmed_pca$sdev^2) / sum(knight_CSF_confirmed_pca$sdev^2)
knight_CSF_confirmed_pca_var_df <- data.frame(component = 1:length(knight_CSF_confirmed_pca_var), 
                                              variance = knight_CSF_confirmed_pca_var)
ggplot(knight_CSF_confirmed_pca_var_df, aes(x = component, y = variance)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 14, color = "red") + 
  labs(x = "Principal Component",
       y = "Proportion of Variance Explained",
       title = "Scree Plot (CSF)") +
  theme_minimal()
knight_plasma_CA_CO <- scale(knight_plasma_prots[index_CA_CO, -c(1:3)], 
      knight_plasma_confirmed_pca$center, knight_plasma_confirmed_pca$scale) %*% 
  knight_plasma_confirmed_pca$rotation 



knight_CSF_confirmed_tsne <- Rtsne(knight_CSF_confirmed[-c(1:3)], check_duplicates = FALSE)
knight_CSF_confirmed_tsne_df <- data.frame(
  X = knight_CSF_confirmed_tsne$Y[, 1],
  Y = knight_CSF_confirmed_tsne$Y[, 2],
  Species = knight_patients_confirmed$CSF_finalstatus
)
ggplot(knight_CSF_confirmed_tsne_df, aes(X, Y, color=Species)) + geom_point() + theme_minimal()

knight_plasma_confirmed_tsne <- Rtsne(knight_plasma_confirmed[-c(1:3)], check_duplicates = FALSE)
knight_plasma_confirmed_tsne_df <- data.frame(
  X = knight_plasma_confirmed_tsne$Y[, 1],
  Y = knight_plasma_confirmed_tsne$Y[, 2],
  Species = knight_patients_confirmed$CSF_finalstatus
)
ggplot(knight_plasma_confirmed_tsne_df, aes(X, Y, color=Species)) + geom_point() + theme_minimal()

knight_ratios_confirmed_tsne <- Rtsne(knight_ratios_confirmed, check_duplicates = FALSE)
knight_ratios_confirmed_tsne_df <- data.frame(
  X = knight_ratios_confirmed_tsne$Y[, 1],
  Y = knight_ratios_confirmed_tsne$Y[, 2],
  Species = knight_patients_confirmed$CSF_finalstatus
)
ggplot(knight_ratios_confirmed_tsne_df, aes(X, Y, color=Species)) + geom_point() + theme_minimal()


perplexities <- c(5, 10, 20, 30, 40, 50, 70, 100, 200, 500)
plots <- list()
for(i in 1:length(perplexities)){
  tsne_result <- Rtsne(knight_CSF_confirmed[-c(1:3)], perplexity = perplexities[[i]])
  tsne_df <- data.frame(
    X = tsne_result$Y[, 1],
    Y = tsne_result$Y[, 2],
    Species = knight_patients_confirmed$CSF_finalstatus
  )
  plots[[i]] <- ggplot(tsne_df, aes(X, Y, color=Species)) + 
    geom_point() + theme_minimal()
}





knight_CSF_confirmed_umap <- umap(knight_CSF_confirmed[-c(1:3)], n_components = 3)
knight_CSF_confirmed_umap_out <- cbind(data.frame(knight_CSF_confirmed_umap$layout), knight_patients_confirmed$CSF_finalstatus)
plot_ly(knight_CSF_confirmed_umap_out, x = ~X1, y = ~X2, z = ~X3, color = ~knight_patients_confirmed$CSF_finalstatus, 
        colors = c('#636EFA','#EF553B','#00CC96')) 

knight_plasma_confirmed_umap <- umap(knight_plasma_confirmed[-c(1:3)], n_components = 3)
knight_plasma_confirmed_umap_out <- cbind(data.frame(knight_plasma_confirmed_umap$layout), knight_patients_confirmed$CSF_finalstatus)
plot_ly(knight_plasma_confirmed_umap_out, x = ~X1, y = ~X2, z = ~X3, color = ~knight_patients_confirmed$CSF_finalstatus, 
        colors = c('#636EFA','#EF553B','#00CC96')) 

knight_ratios_confirmed_umap <- umap(knight_ratios_confirmed[-c(1:3)], n_components = 3)
knight_ratios_confirmed_umap_out <- cbind(data.frame(knight_ratios_confirmed_umap$layout), knight_patients_confirmed$CSF_finalstatus)
plot_ly(knight_ratios_confirmed_umap_out, x = ~X1, y = ~X2, z = ~X3, color = ~knight_patients_confirmed$CSF_finalstatus, 
        colors = c('#636EFA','#EF553B','#00CC96')) 
  

knight_confirmed_lmodel <- generate_lmodels(knight_ratios_confirmed, knight_patients_confirmed,
                                            c("Sex", "plasma_drawage", "CSF_finalstatus"))
knight_confirmed_lmsummary <- generate_lmsummary(knight_confirmed_lmodel, 
                                                 colnames(knight_ratios_confirmed), 
                                                 c("Male", "Age", "AD"))
knight_confirmed_volcanodata <- generate_volcanodata(knight_confirmed_lmsummary, 
                                                     c("Male", "Age", "AD"))
generate_volcanoplot(knight_confirmed_volcanodata, c("Male", "Age", "AD"), 3, 15)
knight_ratios_confirmed_int <- int_dataframe(knight_ratios_confirmed)
knight_confirmed_lmodel_int <- generate_lmodels(knight_ratios_confirmed_int, knight_patients_confirmed,
                                            c("Sex", "plasma_drawage", "CSF_finalstatus"))
knight_confirmed_lmsummary_int <- generate_lmsummary(knight_confirmed_lmodel_int, 
                                                 colnames(knight_ratios_confirmed_int), 
                                                 c("Male", "Age", "AD"))
knight_confirmed_volcanodata_int <- generate_volcanodata(knight_confirmed_lmsummary_int, 
                                                     c("Male", "Age", "AD"))
generate_volcanoplot(knight_confirmed_volcanodata_int, c("Male", "Age", "AD"), 3, 20)

top_AD_genes <- (knight_confirmed_volcanodata[[3]])[order(-knight_confirmed_volcanodata[[3]]$qval), ]
top_AD_upregulated_genes <- top_AD_genes[top_AD_genes$diffexpressed != "none", ]$Protein

knight_ratios_conf_AD_mean <- colMeans(knight_ratios[knight_confirmed_AD, top_AD_upregulated_genes])
knight_ratios_conf_CO_mean <- colMeans(knight_ratios[knight_confirmed_CO, top_AD_upregulated_genes])

lower_CA_CO <- rowSums((knight_ratios[index_CA_CO, top_AD_upregulated_genes] < knight_ratios_conf_CO_mean) / length(top_AD_upregulated_genes))
between_CA_CO <- rowSums((knight_ratios[index_CA_CO, top_AD_upregulated_genes] > knight_ratios_conf_CO_mean & 
                            knight_ratios[index_CA_CO, top_AD_upregulated_genes] < knight_ratios_conf_AD_mean) / length(top_AD_upregulated_genes))
higher_CA_CO <- rowSums((knight_ratios[index_CA_CO, top_AD_upregulated_genes] > knight_ratios_conf_AD_mean) / length(top_AD_upregulated_genes))
data.frame(lower = lower_CA_CO, between = between_CA_CO, higher = higher_CA_CO)


comparison_CA_CO <- knight_ratios[index_CA_CO, top_AD_upregulated_genes] - knight_ratios_conf_CO_mean
rowSums((comparison_CA_CO > 0) / length(top_AD_upregulated_genes)) # proportion of upregulated genes

comparison_CO_AD <- knight_ratios[index_CO_AD, top_AD_upregulated_genes] - knight_ratios_conf_CO_mean
rowSums((comparison_CO_AD > 0) / length(top_AD_upregulated_genes))

index_CA_CO <- knight_patients$CSF_drawstatus == "CA" & knight_patients$CSF_finalstatus == "CO"
patients_CA_CO <- knight_patients[index_CA_CO, ]
index_CO_AD <- knight_patients$CSF_drawstatus == "CO" & knight_patients$CSF_finalstatus == "AD"
patients_CO_AD <- knight_patients[index_CO_AD, ]
index_CO_D <- knight_CSF_diagnosis$at_drawdate == "CO" & knight_CSF_diagnosis$final == "Non-AD dementia"
patients_CO_D <- knight_patients[index_CO_D, ]
index_AD_CO <- knight_plasma_diagnosis$at_drawdate == "AD" & knight_plasma_diagnosis$final == "CO" 
patients_AD_CO <- knight_patients[index_AD_CO, ]


index_CSF_anomaly <- (knight_patients$CSF_drawstatus == "CA" & knight_patients$CSF_finalstatus == "CO" ) |
  (knight_patients$CSF_drawstatus == "CO" & knight_patients$CSF_finalstatus != "CO" )
index_CSF_anomaly <- if_else(is.na(index_CSF_anomaly), FALSE, index_CSF_anomaly)
index_plasma_anomaly <- (knight_patients$plasma_drawstatus == "AD" & knight_patients$plasma_finalstatus != "AD" ) |
  (knight_patients$plasma_drawstatus == "CO" & knight_patients$plasma_finalstatus != "CO" )
index_plasma_anomaly <- if_else(is.na(index_plasma_anomaly), FALSE, index_plasma_anomaly)
index_disagree <- knight_patients$CSF_finalstatus != knight_patients$plasma_finalstatus
index_anomaly <- index_CSF_anomaly | index_plasma_anomaly | index_disagree

knight_patients_anomalous <- knight_patients[index_anomaly, ]
write.csv(knight_patients_anomalous, '~/Developer/proteomics_project/generated_data/knight_anomalous.csv', fileEncoding = "UTF-8", 
          row.names = FALSE)
#early disease = all of knight_patients_anomalous has at least one non-CO entry. This means if we order the drawdates
#and final dates in chronological order, and it goes from CO to AD progression, they would be early AD symptoms. 
# Function to check the sequence
knight_patients_early <- knight_patients_anomalous[
  (knight_patients_anomalous$CSF_finalstatus == "AD" | knight_patients_anomalous$CSF_finalstatus == "CO") 
  & (knight_patients_anomalous$plasma_drawstatus == "AD" | knight_patients_anomalous$plasma_drawstatus == "CO") 
  & (knight_patients_anomalous$plasma_finalstatus == "AD" | knight_patients_anomalous$plasma_finalstatus == "CO"), ]
knight_patients_early <- knight_patients_early[c(1, 3, 4, 5, 6, 7, 9, 12, 13, 14, 15, 16, 17), ]
write.csv(knight_patients_early, '~/Developer/proteomics_project/generated_data/knight_early.csv', fileEncoding = "UTF-8", 
          row.names = FALSE)






check_sequence <- function(s) {
  ad_idx <- which(s == 'AD')[1]
  # If 'AD' was never observed
  if (is.na(ad_idx)) {
    return(FALSE)
  }
  # Split into before and after the first 'AD'
  before_ad <- s[1:(ad_idx-1)]
  after_ad <- s[ad_idx:length(s)]
  
  # Return true if all before are 'CO' and all after (including the AD index) are 'AD'
  return(all(before_ad == 'CO', na.rm = TRUE) & all(after_ad == 'AD', na.rm = TRUE))
}

# Check each column and combine results
df$Result <- mapply(check_sequence, knight_patients_anomalous$X, knight_patients_anomalous$Y, knight_patients_anomalous$Z, knight_patients_anomalous$W) %>% Reduce(`&`, .)

print(df)
index_early_disease <- 

## -----------
stanford_diagnosis <- stanford_patients %>% 
  mutate(after = if_else(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus)) %>% 
  mutate(before = if_else(Plasma_drawdate < CSF_drawdate, Plasma_drawstatus, CSF_drawstatus)) %>% 
  select(before, after) %>% 
  mutate(before = if_else(before == "PDMCI", "MCI-PD", before)) %>%
  mutate(after = if_else(after == "PDMCI", "MCI-PD", after)) %>% 
  mutate(before = if_else(before == "ADMCI", "MCI-AD", before)) %>%
  mutate(after = if_else(after == "ADMCI", "MCI-AD", after))

knight_diagnosis <- knight_patients %>% 
  mutate(after = if_else(plasma_drawdate < CSF_drawdate, CSF_drawstatus, plasma_drawstatus)) %>% 
  mutate(before = if_else(plasma_drawdate < CSF_drawdate, plasma_drawstatus, CSF_drawstatus)) %>% 
  select(before, after) 

knight_freq_table <- as.data.frame(table(knight_diagnosis))
ggplot(data = knight_freq_table,
       aes(axis1 = before, axis2 = after, y = Freq)) +
  geom_alluvium(aes(fill = before)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after)) +
  theme_minimal() +
  labs(title = "Transitions from 'before' to 'after'",
       x = NULL, y = "Count")

stanford_freq_table <- as.data.frame(table(stanford_diagnosis))
ggplot(data = stanford_freq_table,
       aes(axis1 = before, axis2 = after, y = Freq)) +
  geom_alluvium(aes(fill = before)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after)) +
  theme_minimal() +
  labs(title = "Transitions from 'before' to 'after'",
       x = NULL, y = "Count")

knight_CSF_freq_table <- as.data.frame(table(knight_CSF_diagnosis))
ggplot(data = knight_CSF_freq_table,
       aes(axis1 = at_drawdate, axis2 = final, y = Freq)) +
  geom_alluvium(aes(fill = at_drawdate)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = final)) +
  theme_minimal() +
  labs(title = "CSF Transitions from 'at draw' to 'final'",
       x = NULL, y = "Count")

knight_plasma_freq_table <- as.data.frame(table(knight_plasma_diagnosis))
ggplot(data = knight_plasma_freq_table,
       aes(axis1 = at_drawdate, axis2 = final, y = Freq)) +
  geom_alluvium(aes(fill = at_drawdate)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = final)) +
  theme_minimal() +
  labs(title = "Plasma Transitions from 'at draw' to 'final'",
       x = NULL, y = "Count")


# Weird entries
# CSF: 5 CA -> CO, 6 CO -> AD, 7 CO -> Non-AD Dementia
# There's also CA -> DLB, CO -> FTD, etc, but assuming these aren't significant. 
# (uninterested, small #)
# Plasma: 13 CO -> AD, 5 AD -> CO, 6 CO -> Non-AD Dementia, 9 AD -> Non-AD Dementia
# remaining are small or uninterested. 
## CSF: CA -> CO 


# The weird ones --> are they trhe same ones for CSF and plasma? 
# T-SNE or PCA or DEA to see where we should put them. 
