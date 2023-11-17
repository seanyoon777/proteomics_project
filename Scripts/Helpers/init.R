#setting foundations
load_lib <- function(packages, repos = "http://cran.us.r-project.org") {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = repos)
    }
    library(package, character.only = TRUE)
  }
}
#BiocManager::install("preprocessCore")

load_lib(c("sva", "tidyverse", "caret", "data.table", "igraph", "cluster", "purrr", 
           "corrplot", "dplyr", "ggplot2", "ggraph", "circlize", 
           "ggrepel", "openxlsx", "stringr", "tidyr", "WGCNA", "gridExtra"))
load_lib(c("corrr", "gprofiler2", "ComplexHeatmap"))
