source("Scripts/Helpers/init.R")

getProteinNames <- function(proteins) {
  return(str_extract(proteins, "\\w+(?=\\.)"))
}

enrichVolcanodata <- function(data, db) {
  enriched <- enrichr(str_extract(data$Protein, "\\w+(?=\\.)"), db)
  return(enriched)
}
