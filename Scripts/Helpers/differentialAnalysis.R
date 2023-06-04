
generate_lmodels <- function(protdata, patientdata, colnames) {
  lmodels <- list()
  all_prots <- ncol(protdata)
  for(i in 1:ncol(protdata)) {
    new_patientdata <- patientdata[colnames] %>% mutate(protein = protdata[, i])
    lmodels[[i]] <- summary(lm(protein ~ ., data = new_patientdata))
  }
  return(lmodels)
}


generate_lmsummary <- function(lmodels, all_prots, xVars) {
  lmsummary <- data.frame()
  for (i in 1:length(lmodels)) {
    tidyresult <- data.frame(lmodels[[i]]$coefficients) %>% 
      select(Estimate, Std..Error, Pr...t..)
    colnames(tidyresult) <- c("Coefficient", "StdError", "Pval")
    tidyresult <- tidyresult[2:nrow(tidyresult), ] %>% 
      mutate(xVar = xVars, 
             Protein = all_prots[i]) %>% 
      dplyr::select(Protein, xVar, Coefficient, StdError, Pval)
    rownames(tidyresult) <- c()
    lmsummary <- rbind(lmsummary, tidyresult)
  }
  return(lmsummary)
}


generate_volcanodata <- function(data, xVars) {
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
                      size = 3, color = "black")
  }
  
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = 2)
}

getTopProteins <- function(volcanodata, regulated) {
  if (regulated == "up"){
    temp <- volcanodata[volcanodata$diffexpressed == "Upregulated", ] %>% 
      arrange(-qval)
  } 
  if (regulated == "down") {
    temp <- volcanodata[volcanodata$diffexpressed == "Downregulated", ] %>% 
      arrange(-qval)
  }
  return(temp)
}
