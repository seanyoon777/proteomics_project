source("Scripts/Helpers/init.R")

getProteinNames <- function(proteins) {
  return(str_extract(proteins, "\\w+(?=\\.)"))
}

enrichVolcanodata <- function(data, db) {
  enriched <- enrichr(str_extract(data$Protein, "\\w+(?=\\.)"), db)
  return(enriched)
}

enrichClusts <- function(clustdata) {
  enriched <- list()
  for(i in 1:length(unique(clustdata$cluster))) {
    clustproteins <- clustdata[clustdata$cluster == i, ]$protein
    enriched[[i]] <- enrichr(str_extract(clustproteins, "\\w+(?=\\.)"), "ChEA_2022")$ChEA_2022
  }
  return(enriched)
}
