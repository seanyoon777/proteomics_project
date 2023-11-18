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

all_patientdata <- get_biodata("combined/patientdata.csv")
all_patientdata_AD <- get_biodata("combined/patientdata_AD.csv")
all_patientdata_CO <- get_biodata("combined/patientdata_CO.csv")

all_plasma <- get_biodata("combined/plasma_raw_lodfiltered.csv")
all_plasma_AD <- get_biodata("combined/plasma_raw_AD_lodfiltered.csv")
all_plasma_CO <- get_biodata("combined/plasma_raw_CO_lodfiltered.csv")

all_CSF <- get_biodata("combined/CSF_raw_lodfiltered.csv")
all_CSF_AD <- get_biodata("combined/CSF_raw_AD_lodfiltered.csv")
all_CSF_CO <- get_biodata("combined/CSF_raw_CO_lodfiltered.csv")

rank_transform <- function(df) {
  temp <- as.data.frame(scale(df))
  temp <- as.data.frame(t(apply(temp, 1, rank)))
  return(temp)
}

# Experiment 3: ComBat
batch <- factor(ifelse(all_patientdata$batch_effect == 1, "batch1", "batch2"))
all_CSF <- as.data.frame(t(ComBat(dat = t(all_CSF), batch = batch, mod = NULL)))
all_plasma <- as.data.frame(t(ComBat(dat = t(all_plasma), batch = batch, mod = NULL)))
all_CSF_AD <- all_CSF[all_patientdata$final_status == "AD", ]
all_CSF_CO <- all_CSF[all_patientdata$final_status == "CO", ]
all_plasma_AD <- all_plasma[all_patientdata$final_status == "AD", ]
all_plasma_CO <- all_plasma[all_patientdata$final_status == "CO", ]

# Experiment 1: Raw values  
all_CSF_plasma <- all_CSF / all_plasma
all_CSF_plasma_AD <- all_CSF_AD / all_plasma_AD
all_CSF_plasma_CO <- all_CSF_CO / all_plasma_CO

# Experiment 2 default: z-score, rank norm
all_plasma <- scale(all_plasma)
all_CSF <- scale(all_CSF)

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
dea(all_CSF_plasma, all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 30, reverse_bools)
dedata_ratio_all <- dea_without_plot(all_CSF_plasma, all_patientdata, xVars, xLabs_lmsummary)
printTop3(dedata_ratio_all, "GO_Biological_Process_2023", cutoff = 0)


# AD/CO Separated Volcanoplots
xVars <- c("Sex", "avg_drawage", "drawdate_diff", "storage_days", "batch_effect")
xLabs_lmsummary <- c("Male", "Age", "Drawdate difference", "Storage days", "Study bias")
xLabs_plot <- c("Male", "Age")
reverse_bools <- c(FALSE, FALSE)
dedata_ratio_AD <- dea(log10(all_CSF_plasma_AD), all_patientdata_AD, xVars, xLabs_lmsummary, xLabs_plot, 2, 16, reverse_bools)
dedata_ratio_CO <- dea(log10(all_CSF_plasma_CO), all_patientdata_CO, xVars, xLabs_lmsummary, xLabs_plot, 2, 16, reverse_bools)

# CSF & Plasma Volcano plots 
xVars <- c("Sex", "avg_drawage", "final_status", "drawdate_diff", "storage_days", "batch_effect")
xLabs_lmsummary <- c("Male", "Age", "AD", "drawdate_diff", "storage_days", "batch_effect")
xLabs_plot <- c("Male", "Age", "AD")
reverse_bools <- c(FALSE, FALSE, TRUE)
dedata_CSF_all <- dea(log10(all_CSF), all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 16, reverse_bools)
dedata_plasma_all <- dea(log10(all_plasma), all_patientdata, xVars, xLabs_lmsummary, xLabs_plot, 3, 16, reverse_bools)

# Robustness Analysis: Volcanoplots for 20 & 40 day interval drawdate bins
drawdate_bins_20 <- c(0, 1, 20, 40, 60, 80, 100, 120)
dedata_by_20 <- dedata_by_bin(log10(all_CSF_plasma), all_patientdata, drawdate_bins_20)
dea_by_bin(dedata_by_20, drawdate_bins_20, 1, ncol = 4)
pairwise_scatter(dedata_by_20$dedata, drawdate_bins_20, 2, "log2fc")

drawdate_bins_30 <- c(0, 30, 60, 90, 120)
dedata_by_30 <- dedata_by_bin(all_CSF_plasma, all_patientdata, drawdate_bins_30)
pairwise_scatter(dedata_by_30$dedata, drawdate_bins_30, 3, "log2fc")
dea_by_bin(dedata_by_30, drawdate_bins_30, 1, ncol = 2)


drawdate_bins_40 <- c(0, 1, 40, 80, 120)
dedata_by_40 <- dedata_by_bin(log10(all_CSF_plasma), all_patientdata, drawdate_bins_40)
dea_by_bin(dedata_by_40, drawdate_bins_40, 1, ncol = 2)

pairwise_scatter(dedata_by_40$dedata, drawdate_bins_40, 2, "qval")

drawdate_bins_new <- c(0, 30, 55, 75, 120)
dedata_by_new <- dedata_by_bin(log10(all_CSF_plasma), all_patientdata, drawdate_bins_new)
dea_by_bin(dedata_by_new, drawdate_bins_new, 1, ncol = 2)
pairwise_scatter(dedata_by_new$dedata, drawdate_bins_new, 1, "log2fc")

# Robustness Analysis: log2fc vs log2fc for 20 & 40 day interval drawdate bins


# Robustness Analysis: p value vs p value for 20 & 40 day interval drawdate bins




