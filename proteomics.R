#setting foundations

load_lib("tidyverse")
load_lib("caret")
load_lib("data.table")
load_lib("corrr")
load_lib("igraph")
load_lib("cluster")
load_lib("factoextra")
load_lib("purrr")
load_lib("corrplot")
load_lib("dplyr")
load_lib("ggplot2")
load_lib("pls")
load_lib("ggraph")

loess_models <- function(plasma_prot_zscore_age, loess_span = 0.75) {
  plasma_long <- plasma_prot_zscore_age %>%
    gather(key = "protein", value = "zscore", -Age)
  protein_names <- names(plasma_prot_zscore_age)[-1]
  loess_curves <- vector("list", length(protein_names))
  
  for (i in 1:length(protein_names)) {
    x <- plasma_prot_zscore_age$Age
    y <- plasma_prot_zscore_age[, i+1]
    loess_curves[[i]] <- loess(y ~ x)  # You can adjust the span parameter as needed
  }
  return(loess_curves)
}

loess_data <- function(plasma_prot_zscore_age, loess_curves) {
  protein_names <- names(plasma_prot_zscore_age)[-1]
  age_series <- seq(45, 95, by = 1)
  plasma_pred_zscore <- data.frame(Age = age_series)
  
  for (i in 1:length(protein_names)) {
    x <- plasma_prot_zscore_age$Age
    y <- plasma_prot_zscore_age[, i+1]
    pred <- predict(loess_curves[[i]], age_series)
    col_name <- protein_names[i]
    plasma_pred_zscore[[col_name]] <- pred
  }
  
  plasma_pred_zscore_long <- plasma_pred_zscore %>%
    gather(key = "protein", value = "zscore", -Age)
  return(plasma_pred_zscore_long)
}

loess_heatmap <- function(plasma_pred_zscore_long) {
  ggplot(plasma_pred_zscore_long, aes(x = Age, y = protein, fill = zscore)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "yellow") + 
    ggtitle("Protein Trajectories (CSF)") +
    xlab("Age (years)") +
    ylab("4,987 CSF Proteins") +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
}


#cluster by trend 




