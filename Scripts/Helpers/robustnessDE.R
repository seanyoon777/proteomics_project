source("Scripts/Helpers/differentialAnalysis.R")

dea <- function(protdata, patientdata, xVars, xLabs_lmsummary, xLabs_plot, 
                ncol_plot, prot_num, reverse_bools) {
  lmodels <- generate_lmodels(protdata, patientdata, xVars)
  lsummary <- generate_lmsummary(lmodels, colnames(protdata), xLabs_lmsummary)
  volcanodata <- generate_volcanodata(lsummary, xLabs_plot)
  generate_volcanoplot(volcanodata, xLabs_plot, ncol_plot, prot_num, reverse_bools)
  return(volcanodata)
}

dea_without_plot <- function(protdata, patientdata, xVars, xLabs_lmsummary) {
  lmodels <- generate_lmodels(protdata, patientdata, xVars)
  lsummary <- generate_lmsummary(lmodels, colnames(protdata), xLabs_lmsummary)
  volcanodata <- generate_volcanodata(lsummary, xLabs_plot)
  return(volcanodata)
}

dedata_by_bin <- function(protdata, patientdata, drawdate_bins) {
  volcanodata <- list()
  nsamples <- list()
  num <- length(drawdate_bins) - 1
  
  for(i in 1:num) {
    idx <- patientdata$drawdate_diff < drawdate_bins[i + 1] & patientdata$drawdate_diff >= drawdate_bins[i]
    new_patientdata <- all_patientdata[idx, ]
    new_protdata <- protdata[idx, ]
    
    if(sum(new_patientdata$batch_effect == 1) * sum(new_patientdata$batch_effect == 0) == 0) {
      lmodels <- generate_lmodels(new_protdata, new_patientdata, c("Sex", "avg_drawage", "final_status"))
      lmsummary <- generate_lmsummary(lmodels, colnames(protdata), c("Male", "Age", "AD"))
    } else {
      lmodels <- generate_lmodels(new_protdata, new_patientdata, c("Sex", "avg_drawage", "final_status", "batch_effect"))
      lmsummary <- generate_lmsummary(lmodels, colnames(protdata), c("Male", "Age", "AD", "batch_effect"))
    }
    volcanodata[[i]] <- generate_volcanodata(lmsummary, c("Male", "Age", "AD"))
    nsamples[[i]] <- nrow(new_patientdata)
  }
  res <- list(dedata = volcanodata, nsamples = nsamples)
  return(res)
}

dea_by_bin <- function(res, drawdate_bins, factor_num) {
  titles <- c('Sex', 'Aging', 'Alzheimers')
  reverse_bools <- c(c(FALSE), c(FALSE), c(TRUE))
  num <- length(drawdate_bins) - 1
  volcanoplot <- list()
  for(i in 1:num) {
    text_label <- paste0(drawdate_bins[i], "-", drawdate_bins[i + 1], " days, ", res$nsamples[[i]] , " samples")
    volcanoplot[[i]] <- generate_volcanoplot(res$dedata[[i]], text_label, 1, prot_num, reverse_bools[factor_num])
  }
  
  title <- paste0(titles[factor_num], " Proteome by Drawdate Difference")
  title_grob <- textGrob(title, gp=gpar(fontsize=20, fontface="bold"))
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = ncol, top = title_grob)
}


