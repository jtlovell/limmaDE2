#' @title Find modules that have high connectivity
#'
#' @description
#' \code{modulesByConnectivity} Determine the strongest modules in a network
#' using the mean absolute value of signedKME.
#'
#' @param net WGCNA generated network. Usually from WGCNA::blockwiseModules
#' @param datExpr The expression dataset, transposed so that genes are columns
#' and individuals are rows.
#' @param mean.kme.threshold Numeric, the thresholds to test.
#' Can be of length >1.
#' @details More here soon.

#' @return A vector of modules for each kme threshold specified.
#' @examples
#' \dontrun{
#' library(WGCNA)
#' library(igraph)
#' data(kidney) #' from simseq
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-with(kidney,
#'            data.frame(id = paste(replic, treatment, sep = "_"),
#'                       rep=replic,
#'                       Treatment=ifelse(treatment == "Tumor","tumor","cntr"),
#'                       stringsAsFactors=F))
#' colnames(counts)<-info$id
#' stats <- pipeLIMMA(counts = counts,
#'                    info = info,
#'                    block = NULL,
#'                    formula = "~ Treatment")
#'
#' datExpr.1=t(stats$voom$E)
#' pow=6
#' net.1 = blockwiseModules(datExpr.1, power = pow,
#'                          maxBlockSize = 10000, deepSplit = 2,
#'                          minModuleSize = 10,
#'                          saveTOMs = FALSE,
#'                          verbose = F)
#' modulesByConnectivity(net = net.1, datExpr = datExpr.1)
#' }
#' @export
modulesByConnectivity<-function(datExpr,net, mean.kme.threshold = seq(0,1,.1)){
  datKME=signedKME(datExpr, datME = net$MEs, outputColumnName = "")
  mods = net[[1]]

  mean.kme<-sapply(unique(mods), function(x){
    mean(abs(datKME[mods==x,x]))
  })
  good.mods<-lapply(mean.kme.threshold, function(x){
    names(mean.kme)[mean.kme>x]
  })
  names(good.mods)<-paste0("kme_",mean.kme.threshold)
  return(good.mods)
}
