
generate_lmodels <- function(protdata, patientdata, colnames) {
  lmodels <- list()
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
             log2fc = Coefficient, se = StdError) %>% 
      select(Protein, Pval, padj, qval, se, log2fc) %>% 
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
      select(Protein, Pval, padj, qval, se, log2fc) %>% 
      mutate(diffexpressed = ifelse(log2fc > 0 & qval >= 1,
                                    yes = "Upregulated", 
                                    no = ifelse(log2fc < 0 & qval >= 1,
                                                yes = "Downregulated", no="none")))
  }
  return(volcano_data)
}


generate_volcanoplot <- function(volcanodata, xVars, ncols, prot_nums) {
  volcanoplot <- list()
  for(i in 1:length(xVars)) {
    
    volcanodata_temp <- volcanodata[[i]]
    top_volcanodata <- volcanodata_temp[volcanodata_temp$diffexpressed != "none", ]
    top_genes <- head((volcanodata[[i]])[order(-volcanodata[[i]]$qval), ], prot_nums)
    max_qval <- max(abs(volcanodata_temp$qval), na.rm = TRUE)
    
    volcanodata_temp$point_size <- abs(volcanodata_temp$qval)^2 / max_qval
    if (sum(volcanodata_temp$diffexpressed != "none") == 0) {
      volcanodata_temp$point_size <- 0.1
      scale_size <- c(0.5, 1.5)
    } else {
      scale_size <- c(1, 6)
    }
    
    if (xVars[i] == "AD") {
      volcanoplot[[i]] <- ggplot(volcanodata_temp, aes(x = -log2fc, y = qval, color = factor(diffexpressed))) + 
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
        scale_color_manual(values = c("Downregulated" = "indianred1",
                                      "Upregulated" = "royalblue1",
                                      "none" = "#2c2c2c")) +
        annotate("text", x = min((volcanodata[[i]])$log2fc)*1.1, y = max((volcanodata[[i]])$qval), 
                 label = xVars[i], hjust = 0, size = 4) +
        geom_text_repel(data = top_genes, aes(label = Protein, y = qval, x = -log2fc),
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
    } else {
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
        annotate("text", x = min((volcanodata[[i]])$log2fc)*1.1, y = max((volcanodata[[i]])$qval), 
                 label = xVars[i], hjust = 0, size = 4) +
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
  
  }
  
  gridExtra::grid.arrange(grobs = volcanoplot, ncol = ncols)
}

inverse_normal_transform <- function(x) {
  ranks <- rank(x)
  percentiles <- ranks / (length(x) + 1)
  transformed <- qnorm(percentiles)
  return(transformed)
}

int_dataframe <- function(data) {
  df <- data.frame(apply(data, 2, inverse_normal_transform))
  return(df)
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

correlationplot <- function(volcanodata, varindex) {
  interval <- 120 / length(volcanodata) 
  titles <- c("Sex", "Aging", "Alzheimers")
  layout_matrix <- matrix(NA, nrow = length(volcanodata), ncol = length(volcanodata))
  cnt <- 1
  for (i in 1:(length(volcanodata))) {
    for (j in i:(length(volcanodata))) {
      layout_matrix[i, j] <- cnt
      cnt <- cnt + 1
    }
  }
  vp <- viewport(layout = grid.layout(nrow = length(volcanodata) + 1, 
                                      ncol = length(volcanodata) + 1))
  
  pushViewport(vp)
  cnt <- 2
  bins <- c(0, interval * 0:(length(volcanodata) - 1))
  
  for (i in 2:(length(volcanodata) + 1)) {
    print(textGrob(label = paste0(bins[i - 1], '-', bins[i], ' days'), gp = gpar(fontsize = 12)), 
          vp = viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  for (i in 2:(length(volcanodata) + 1)) {
    print(textGrob(label = paste0(bins[i - 1], '-', bins[i], ' days'), gp = gpar(fontsize = 12)), 
          vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  
  
  volcanoplots <- list()
  for (i in 1:(length(volcanodata) - 1)) {
    for (j in 2:length(volcanodata)) {
      if(j > i) {
        data <- data.frame(x = volcanodata[[i]][[varindex]]$log2fc, y = volcanodata[[j]][[varindex]]$log2fc)
        axismax <- max(max(data$x), max(data$y))
        axismin <- min(min(data$x), min(data$y))
        volcanoplots[[cnt]] <- ggplotGrob(ggplot(data, aes(x = x, y = y)) +
          geom_point(color = "steelblue", alpha = 0.3, size = 2) + # Set color and alpha
          theme_minimal() +
          labs(x = "X Values", y = "Y Values") + 
          xlim(c(axismin, axismax)) + 
          ylim(c(axismin, axismax)) + 
          annotate("text", x = mean(data$x), y = mean(data$y), label = paste0("i=", i, ", j=", j)) +
          geom_smooth(method = "lm", aes(color = 'red')) +
          theme(strip.text = element_blank(), 
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(), 
                axis.text = element_text(size = 12), axis.title = element_blank(), 
                legend.position = "none") 
        )
        print(volcanoplots[[cnt]], vp = viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
        cnt <- cnt + 1
      }
    }
  }
  popViewport()

  # Add row and column labels

}


grid.text(titles[varindex], gp = gpar(fontsize = 20, fontface = "bold"), y = 1)

row_labels = c('f', 'f', 'e', 'q', 'g', 'a', 'q')
col_labels = c('f', 'f', 'e', 'q', 'g', 'a', 'q')
row_table <- tableGrob(t(as.matrix(row_labels)), theme = ttheme_minimal())
col_table <- tableGrob(as.matrix(col_labels), theme = ttheme_minimal())
grobs <- c(list(row_table, col_table), volcanoplots)
layout_matrix[1, ] <- 1
layout_matrix[, 1] <- 2


title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
title_grob <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"))
grid.arrange(grobs = do.call(grid::gList, grobs), layout_matrix = layout_matrix, top = title_grob)
  title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
  title_grob <- textGrob(title, gp=gpar(fontsize=20, fontface="bold"))
  grid.draw(title_grob)



row_table <- tableGrob(t(as.matrix(row_labels)), theme = ttheme_minimal())
col_table <- tableGrob(as.matrix(col_labels), theme = ttheme_minimal())

grobs <- c(list(row_table, col_table), volcanoplots)
layout_matrix[1, ] <- 1
layout_matrix[, 1] <- 2

title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
title_grob <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"))

grid.arrange(grobs = do.call(grid::gList, grobs), layout_matrix = layout_matrix, top = title_grob)




grobs <- lapply(volcanoplots, ggplotGrob)
title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
title_grob <- textGrob(title, gp=gpar(fontsize=20, fontface="bold"))
grid.arrange(grobs = do.call(grid::gList, grobs), layout_matrix = layout_matrix, top = title_grob)




correlationplot <- function(volcanodata, varindex, interval) {
  grid.newpage()
  bins <- c(0, interval * 0:(length(volcanodata) - 1))

  data <- data.frame()
  cnt <- 1
  for (i in 1:(length(volcanodata) - 1)) {
    for (j in (i + 1):length(volcanodata)) {
      temp <- data.frame(x = volcanodata[[i]][[varindex]]$log2fc, y = volcanodata[[j]][[varindex]]$log2fc)
      temp$row <- i
      temp$col <- j
      data <- rbind(data, temp)
      cnt <- cnt + 1
    }
  }
  
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(color = "steelblue", alpha = 0.3, size = 2) +
    geom_smooth(method = "lm", aes(color = 'red')) +
    facet_grid(row ~ col, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_blank(),
          legend.position = "none")
  
  g <- ggplotGrob(p)
  g <- gtable::gtable_add_rows(g, heights = unit(1, "lines"), 0)
  g <- gtable::gtable_add_cols(g, widths = unit(2, "lines"), 0)
  rownames <- c("Row Labels", as.character(1:(length(volcanodata) - 1)))
  colnames <- c("Col Labels", as.character(2:length(volcanodata)))
  
  grid.draw(g)
}

g <- gtable::gtable_add_grob(g, grobs = textGrob(bins, vjust = 1, hjust = 0), t = 1:length(bins), l = 1)
g <- gtable::gtable_add_grob(g, grobs = textGrob(bins, vjust = 0, hjust = 1, rot = 90), t = 1, l = 1:length(bins))



correlationplot <- function(volcanodata, varindex) {
  titles <- c("Sex", "Aging", "Alzheimers")
  n <- length(volcanodata)
  volcanoplots <- list()
  cnt <- 1
  for (i in 1:n) {
    for (j in i:n) {
      if (j > i) {
        data <- data.frame(x = volcanodata[[i]][[varindex]]$log2fc, y = volcanodata[[j]][[varindex]]$log2fc)
        axismax <- max(max(data$x), max(data$y))
        axismin <- min(min(data$x), min(data$y))
        volcanoplots[[cnt]] <- ggplot(data, aes(x = x, y = y)) +
          geom_point(color = "steelblue", alpha = 0.3, size = 2) +
          theme_minimal() +
          labs(x = "X Values", y = "Y Values") +
          xlim(c(axismin, axismax)) +
          ylim(c(axismin, axismax)) +
          geom_smooth(method = "lm", aes(color = 'red')) +
          theme(strip.text = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.text = element_text(size = 12), legend.position = "none")
        cnt <- cnt + 1
      }
    }
  }
  
  layout_matrix <- matrix(NA, nrow = n, ncol = n)
  cnt <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      layout_matrix[i, j] <- cnt
      cnt <- cnt + 1
    }
  }
  title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
  title_grob <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.arrange(grobs = c(title_grob, lapply(volcanoplots, ggplotGrob)), layout_matrix = rbind(c(0, rep(NA, n - 1)), layout_matrix))
}



correlationplot <- function(volcanodata, varindex, interval) {
  titles <- c("Sex", "Aging", "Alzheimers")
  n <- length(volcanodata) - 1
  bins <- c(0, interval * 0:n)
  volcanoplots <- list()
  cnt <- 1
  for (i in 1:n) {
    for (j in (i + 1):(n + 1)) {
      if (j > i) {
        data <- data.frame(x = rank(ratio_stability_volcanodata_0[[i]][[varindex]]$log2fc 
                                         / ratio_stability_volcanodata_0[[i]][[varindex]]$se), 
                           y = rank(ratio_stability_volcanodata_0[[j]][[varindex]]$log2fc 
                                         / ratio_stability_volcanodata_0[[j]][[varindex]]$se))
        axismax <- max(max(data$x), max(data$y))
        axismin <- min(min(data$x), min(data$y))
        
        lm_model <- lm(y ~ x, data = data)
        r2 <- summary(lm_model)$r.squared
        cor_pearson <- cor(data$x, data$y)
        
        volcanoplots[[cnt]] <- ggplot(data, aes(x = x, y = y)) +
          geom_point(color = "steelblue", alpha = 0.3, size = 2) +
          theme_minimal() +
          xlim(c(axismin, axismax)) +
          ylim(c(axismin, axismax)) +
          geom_smooth(method = "lm", aes(color = 'red'), SE = TRUE, fullrange = TRUE) +
          expand_limits(x = c(axismin, axismax), y = c(axismin, axismax)) +
          annotate("text", x = Inf, y = Inf, vjust = 1.8, hjust = 1, label = paste("R^2: ", round(r2, 3)), size = 3) +
          annotate("text", x = Inf, y = Inf, vjust = 3.1, hjust = 1, label = paste("Pearson: ", round(cor_pearson, 3)), size = 3) +
          theme(strip.text = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.title=element_text(size=14),
                axis.text = element_text(size = 12), legend.position = "none") + 
          labs(x = NULL, y = NULL)
        if (j == n + 1) {
          volcanoplots[[cnt]] <- volcanoplots[[cnt]] + labs(x = paste0(bins[i], '-', bins[i + 1], ' days'))
        } 
        if (i == 1) {
          volcanoplots[[cnt]] <- volcanoplots[[cnt]] + labs(y = paste0(bins[j], '-', bins[j + 1], ' days'))
        }
        cnt <- cnt + 1
      }
    }
  }
  
  cnt <- 1
  
  layout_matrix <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      layout_matrix[j, i] <- cnt
      cnt <- cnt + 1
    }
  }
  
  title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
  title_grob <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"), just = "center")
  col_widths <- rep(1, n)
  col_widths[1] <- 1.1
  row_heights <- rep(1, n)
  row_heights[n] <- 1.1
  
  grobs <- list()
  for (i in 1:length(volcanoplots)) {
    grobs[[i]] <- ggplotGrob(volcanoplots[[i]])
  }
  
  grid.arrange(grobs = grobs, layout_matrix = layout_matrix, widths = col_widths, heights = row_heights)
}
my_plot <- correlationplot(ratio_stability_volcanodata_0, 1, 20) 
correlationplot(ratio_stability_volcanodata_0, 2, 20)
correlationplot(ratio_stability_volcanodata_0, 3, 20)



volcanodata[[1]]


correlationplot_pval <- function(volcanodata, varindex, interval) {
  titles <- c("Sex", "Aging", "Alzheimers")
  n <- length(volcanodata) - 1
  bins <- c(0, interval * 0:n)
  volcanoplots <- list()
  cnt <- 1
  for (i in 1:n) {
    for (j in (i + 1):(n + 1)) {
      if (j > i) {
        data <- data.frame(x = rank(volcanodata[[i]][[varindex]]$Pval), 
                           y = rank(volcanodata[[j]][[varindex]]$Pval))
        axismax <- max(max(data$x), max(data$y))
        axismin <- min(min(data$x), min(data$y))
        
        lm_model <- lm(y ~ x, data = data)
        r2 <- summary(lm_model)$r.squared
        cor_pearson <- cor(data$x, data$y)
        
        volcanoplots[[cnt]] <- ggplot(data, aes(x = x, y = y)) +
          geom_point(color = "steelblue", alpha = 0.3, size = 2) +
          theme_minimal() +
          xlim(c(axismin, axismax)) +
          ylim(c(axismin, axismax)) +
          geom_smooth(method = "lm", aes(color = 'red'), SE = TRUE, fullrange = TRUE) +
          expand_limits(x = c(axismin, axismax), y = c(axismin, axismax)) +
          annotate("text", x = Inf, y = Inf, vjust = 1.8, hjust = 1, label = paste("R^2: ", round(r2, 3)), size = 3) +
          annotate("text", x = Inf, y = Inf, vjust = 3.1, hjust = 1, label = paste("Pearson: ", round(cor_pearson, 3)), size = 3) +
          theme(strip.text = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                axis.title=element_text(size=14),
                axis.text = element_text(size = 12), legend.position = "none") + 
          labs(x = NULL, y = NULL)
        if (j == n + 1) {
          volcanoplots[[cnt]] <- volcanoplots[[cnt]] + labs(x = paste0(bins[i], '-', bins[i + 1], ' days'))
        } 
        if (i == 1) {
          volcanoplots[[cnt]] <- volcanoplots[[cnt]] + labs(y = paste0(bins[j], '-', bins[j + 1], ' days'))
        }
        cnt <- cnt + 1
      }
    }
  }
  
  cnt <- 1
  
  layout_matrix <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      layout_matrix[j, i] <- cnt
      cnt <- cnt + 1
    }
  }
  
  title <- paste0(titles[varindex], " Proteome Correlation (log2fc)")
  title_grob <- textGrob(title, gp = gpar(fontsize = 20, fontface = "bold"), just = "center")
  col_widths <- rep(1, n)
  col_widths[1] <- 1.1
  row_heights <- rep(1, n)
  row_heights[n] <- 1.1
  
  grobs <- list()
  for (i in 1:length(volcanoplots)) {
    grobs[[i]] <- ggplotGrob(volcanoplots[[i]])
  }
  
  grid.arrange(grobs = grobs, layout_matrix = layout_matrix, widths = col_widths, heights = row_heights)
}
correlationplot_pval(ratio_stability_volcanodata_0, 1, 20)

correlationplot(ratio_stability_volcanodata_40_0, 2, 40)


ggsave("plot.png", plot = my_plot, width = 16, height = 14)

grid.arrange(grobs = grobs, layout_matrix = rbind(c(0, rep(NA, n - 1)), layout_matrix), heights = rep(1, n + 1)) # Set uniform sizes


