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

dea_by_bin <- function(res, drawdate_bins, factor_num, ncol, prot_nums = 10) {
  titles <- c('Sex', 'Aging', 'Alzheimers')
  reverse_bools <- c(c(FALSE), c(FALSE), c(TRUE))
  num <- length(drawdate_bins) - 1
  volcanoplot <- list()
  for(i in 1:num) {
    text_label <- paste0(drawdate_bins[i], "-", drawdate_bins[i + 1], " days, ", res$nsamples[[i]], " samples")
    volcanoplot[[i]] <- generate_volcanoplot(res$dedata[[i]], text_label, 1, prot_nums, 
                                             reverse_bools, xVars_choice = factor_num:factor_num)
  }
  
  title <- paste0(titles[factor_num], " Proteome by Drawdate Difference")
  title_grob <- textGrob(title, gp=gpar(fontsize=20, fontface="bold"))
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = ncol, top = title_grob)
}

pairwise_scatter <- function(res, drawdate_bins, factor_num, type) {
  suppressMessages({pairwise_scatter_inner(res, drawdate_bins, factor_num, type)})
}
  
pairwise_scatter_inner <- function(res, drawdate_bins, factor_num, type) {
  plots <- list()
  type_x <- paste0(type, "_x")
  type_y <- paste0(type, "_y")
  num <- length(drawdate_bins) - 1
  for (j in 1:num) {
    for (i in 1:num) {
      if (i != j) {
        merged_data <- merge(res[[i]][[factor_num]], res[[j]][[factor_num]], 
                             by = "Protein", suffixes = c("_x", "_y")) 
        x_values <- merged_data %>% select(all_of(type_x))
        y_values <- merged_data %>% select(all_of(type_y))
        merged_data <- data.frame(Protein = merged_data$Protein, x = x_values[[1]], y = y_values[[1]])
        fit <- lm(y ~ x, data = merged_data)
        r_squared <- summary(fit)$r.squared
        spearman_cor <- cor(merged_data$x, merged_data$y, method = "spearman")
        plots[[paste(i, j)]] <- ggplot(merged_data, aes(x = x, y = y)) +
          geom_point(color = "skyblue", alpha = 0.5) +
          geom_smooth(method = "lm", color = "red", se = TRUE) +
          annotate("text", x = min(merged_data$x), y = max(merged_data$y), 
                   label = sprintf("R² = %.2f\nSpearman = %.2f", r_squared, spearman_cor), 
                   hjust = 0, vjust = 1, size = 3, color = "black") +
          theme_minimal() +
          scale_x_continuous(position = "top") +
          theme(panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.text = element_text(size = 5),
                legend.position = "none") 
      } else {
        values <- res[[i]][[factor_num]] %>% select(all_of(type))
        range <- max(values) - min(values)
        df <- data.frame(Proteins = res[[i]][[factor_num]]$Protein, x = values[[1]])
        plots[[paste(i, j)]] <- ggplot(df, aes(x = x)) +
          geom_histogram(binwidth = range / 100, fill = "lightblue", color = "lightblue") +
          theme_minimal() +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.text = element_text(size = 5),
                axis.title = element_blank(),
                legend.position = "none") + 
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
      
      if (j != 1 & i != j) {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (i != 1 & i != j) {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
      if (j == 1) {  # Top row
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          labs(x = paste0(drawdate_bins[i], "-", drawdate_bins[i + 1], " days")) + 
          theme(axis.title = element_text(size = 3))
      } else {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.title.x = element_blank(), axis.title = element_blank())
      }
      if (i == 1) {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          labs(y = paste0(drawdate_bins[j], "-", drawdate_bins[j + 1], " days"))
      } else {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.title.y = element_blank(), axis.title = element_blank())
      }
    }
  }
  do.call(gridExtra::grid.arrange, c(plots, ncol = num))
}



pairwise_scatter_inner <- function(res, drawdate_bins, factor_num, type) {
  plots <- list()
  type_x <- paste0(type, "_x")
  type_y <- paste0(type, "_y")
  num <- length(drawdate_bins) - 1
  for (j in 1:num) {
    for (i in 1:num) {
      x_title <- paste0(drawdate_bins[i], "-", drawdate_bins[i + 1], " days")
      y_title <- paste0(drawdate_bins[j], "-", drawdate_bins[j + 1], " days")
      if (i != j) {
        merged_data <- merge(res[[i]][[factor_num]], res[[j]][[factor_num]], 
                             by = "Protein", suffixes = c("_x", "_y")) 
        merged_data <- data.frame(Protein = merged_data$Protein, x = merged_data[[type_x]], y = merged_data[[type_y]])
        fit <- lm(y ~ x, data = merged_data)
        r_squared <- summary(fit)$r.squared
        spearman_cor <- cor(merged_data$x, merged_data$y, method = "spearman")
        
        plots[[paste(i, j)]] <- ggplot(merged_data, aes(x = x, y = y)) +
          geom_point(color = "skyblue", alpha = 0.5) +
          geom_smooth(method = "lm", color = "red", se = TRUE) +
          annotate("text", x = min(merged_data$x), y = max(merged_data$y), 
                   label = sprintf("R² = %.2f\nSpearman = %.2f", r_squared, spearman_cor), 
                   hjust = 0, vjust = 1, size = 3, color = "black") +
          theme_minimal() +
          scale_x_continuous(position = "top") +
          theme(panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.text = element_text(size = 7),
                legend.position = "none",
                axis.title.x = element_text(size = 11),  # Ensure axis titles are visible
                axis.title.y = element_text(size = 11)) 
      } else {
        values <- res[[i]][[factor_num]] %>% select(all_of(type))
        range <- max(values) - min(values)
        df <- data.frame(Proteins = res[[i]][[factor_num]]$Protein, x = values[[1]])
        plots[[paste(i, j)]] <- ggplot(df, aes(x = x)) +
          geom_histogram(binwidth = range / 100, fill = "lightblue", color = "lightblue") +
          theme_minimal() +
          scale_x_continuous(position = "top") +
          theme(panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.text = element_text(size = 5),
                legend.position = "none",
                axis.title.x = element_text(size = 11),  # Ensure axis titles are visible
                axis.title.y = element_text(size = 11), 
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      }
      
      if (j == 1) {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + labs(x = x_title) 
      } else {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.title.x = element_blank(), 
                axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (i == 1) {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + labs(y = y_title) 
      } else {
        plots[[paste(i, j)]] <- plots[[paste(i, j)]] + 
          theme(axis.title.y = element_blank(),
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
    }
  }
  do.call(gridExtra::grid.arrange, c(plots, ncol = num))
}



