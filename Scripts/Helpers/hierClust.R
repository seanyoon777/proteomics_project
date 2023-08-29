source("Scripts/Helpers/init.R")

loess_predict <- function(x, Y, x_seq, xlab) {
  models <- vector("list", ncol(Y))
  data_pred <- data.frame(xlab = x_seq)
  for (i in 1:ncol(Y)) {
    models[[i]] <- loess(Y[, i] ~ x)
    data_pred[, i + 1] <- predict(models[[i]], x_seq)
  }
  colnames(data_pred) <- c(xlab, colnames(Y))
  return(data_pred)
}


# Input: hclust object, distance matrix
# Output: silhouette scores of each matrix
# Example --- 
# ratio_CO_z_dist <- dist(t(predratio_CO_z[-1]), method = "euclidean")
# ratio_CO_z_clust <- hclust(ratio_CO_z_dist, method = "complete")
determine_numclust <- function(clust, dist, maxnum = 10) {
  silhouette <- vector()
  for (k in 2:maxnum) {
    cut <- cutree(clust, k)
    vals <- silhouette(cut, dist)
    silhouette[k - 1] <- mean(vals[, 3])
  }
  plot(silhouette)
}

# Returns 
proteinByClust <- function(clust, dist, nclust) {
  clust <- cutree(clust, k = nclust) %>% 
    data.frame(cluster = .)
  clust["protein"] = rownames(clust)
  clust <- clust %>% relocate(cluster, .after = "protein")
  rownames(clust) <- c()
  return(clust)
}

proteinClustData <- function(data, clust) { 
  long <- data %>% gather(key = protein, value = value, -Age) %>% 
    inner_join(clust, by = "protein") %>% 
    arrange(cluster, Age, value)
  return(long)
}

generate_clusterplot <- function(data) {
  data %>% 
    ggplot(aes(x = Age, y = value, group = protein, color = factor(cluster))) +
    geom_line(stat = "smooth", method = "loess", se = FALSE, linewidth = 0.5, alpha = 0.1) +
    scale_color_manual(
      values = c("#FF6F61", "#6B5B95", "#88B04B", "#FFA500", "#92A8D1", "#FF69B4", 
                 "#955251", "#008080", "#DA70D6", "#6A5ACD")) +
    labs(x = "Age (years)", y = "CSF Plasma Protein Ratio (Z-Scored)", color = "Cluster") +
    facet_wrap(~ cluster, ncol = 3, nrow = 4, labeller = labeller(cluster = as.character)) +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.8, alpha = 1, linetype = "solid", color = "white") +
    stat_summary(fun.data = "mean_cl_normal", geom = "line",
                 aes(group = cluster, color = factor(cluster)),
                 linewidth = 0.5, alpha = 0.8) + 
    theme(strip.text = element_blank(), 
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_blank(),
          axis.line = element_blank(), 
          legend.position = "none") +
    labs(title = NULL)
}
