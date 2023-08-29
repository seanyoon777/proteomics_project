source("Scripts/Helpers/init.R")

decide_soft_threshold <- function(data) {
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

cluster_dendrogram <- function(data, soft_threshold, varname, width = 12, height = 9) {
  # Soft thresholding power of 6
  net <- blockwiseModules(data, power = soft_threshold,
                         TOMType = "unsigned", minModuleSize = 20,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "allPlasmaTOM",
                         verbose = 3)
  
  sizeGrWindow(width = width, height = height)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                      main = paste0(varname, " Cluster Dendrogram"))
  return(net)
}

cluster_heatmap_dendrogram <- function(data, net, varname, power) {
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  
  # Define numbers of genes and samples
  nGenes = ncol(data);
  nSamples = nrow(data);
  
  dissTOM = 1-TOMsimilarityFromExpr(data, power = power);
  plotTOM = dissTOM^(power + 1);
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
  # Call the plot function
  sizeGrWindow(9,9)
  TOMplot(plotTOM, geneTree, moduleColors, main = paste0(varname, " Network Heatmap Plot"))
  
  # Recalculate module eigengenes
  MEs = moduleEigengenes(data, moduleColors)$eigengenes
  # Isolate weight from the clinical traits
  weight = as.data.frame(datTraits$weight_g);
  names(weight) = "weight"
}


