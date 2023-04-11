#setting foundations
load_lib <- function(packages, repos = "http://cran.us.r-project.org") {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = repos)
    }
    library(package, character.only = TRUE)
  }
}

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
load_lib("circlize")
load_lib("Cairo")
load_lib("ComplexHeatmap")
load_lib("loessclust")


