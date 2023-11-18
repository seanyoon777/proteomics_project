source("Scripts/Helpers/init.R")

getProteinNames <- function(proteins) {
  return(str_extract(proteins, "\\w+(?=\\.)"))
}

enrichVolcanodata <- function(data, type, db, cutoff = 0) {
  newdata <- data[data$diffexpressed == type & abs(data$log2fc) > cutoff, ]
  enriched <- enrichr(str_extract(newdata$Protein, "\\w+(?=\\.)"), db)
  return(enriched)
}

printTop3 <- function(data, db, cutoff = 0) {
  n <- length(data)
  type <- c("Upregulated", "Downregulated")
  for (i in 1:n) {
    for (j in 1:2) {
      enriched <- enrichVolcanodata(data[[i]], type[j], db, cutoff)
      print(enriched[[1]]$Term[c(1, 2, 3)])
      #out <- paste0("[1] ", enriched[[1]]$Term[1], "\n", "[2] ", enriched[[1]]$Term[2], "\n", "[3]", enriched[[1]]$Term[3])
    }
  }
}

enrichClusts <- function(clustdata) {
  enriched <- list()
  for(i in 1:length(unique(clustdata$cluster))) {
    clustproteins <- clustdata[clustdata$cluster == i, ]$protein
    enriched[[i]] <- enrichr(str_extract(clustproteins, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022
  }
  return(enriched)
}
