#' @title Cannonical correspondence analysis (CCA) of gene expression
#'
#' @description
#' \code{de2cca} A combined multivariate analysis of gene expression that invokes both
#' canonical correspondence analysis of the expression matrix and a hierarchical clustering
#' of the distance matrix.
#' @param counts A count matrix
#' @param info An experimental design matrix
#' @param formula A character string that can be coerced to a formula. Terms must be separated
#' by a whitespace
#' @param fast If the dataset has >10k genes, this indicates that the expression matrix should be
#' subsampled down to 10k. Datasets with >10k genes can take >10 minutes to run
#' @param palette A R palette that is used to color the dendrogram
#'
#' @return a vegan::cca model object
#' @examples
#' library(PhyGenomicsData.PV2016)
#' data(info12); data(counts12)
#' info<-info12
#' counts<-counts12
#' test<-de2cca(info=info, counts=counts[1:1000,], formula="Treatment + MDWP + order")
#' @importFrom  vegan cca
#' @importFrom  dendroextras colour_clusters
#' @export
de2cca<-function(info, counts, fast=TRUE, formula, palette=NULL,...){
  if(is.null(sample) &  nrow(counts)>10000 & fast){
    cat("randomly subsampling genes to 12k, to speed up analysis\n")
    counts<-counts[sample(1:nrow(counts),10000),]
  }

  ind<-complete.cases(info)
  i<-info[ind,]
  c<-counts[,ind]
  grp<-strsplit(formula," ")[[1]][1]
  grp<-i[,grp]
  if(is.null(palette)){
    if(is.numeric(grp)){
      pal<-colorRampPalette(c("red", "white", "blue"))
      pal<-pal(length(unique(grp)))
    }else{
      pal<-rainbow(length(levels(grp)))
    }
  }else{
    if(is.numeric(grp)){
      pal<-colorRampPalette(c("red", "white", "blue"))
      pal<-pal(length(unique(grp)))
    }else{
      pal<-colors
    }
  }
  hc<-as.dendrogram(hclust(dist(t(c))))
  cols<-pal[as.integer(grp)]
  names(cols)<-labels(hc)
  tp<-set_leaf_colours(hc, col=cols, col_to_set = "label")
  plot(tp, main=paste("distance-based hclust dendrogram \n colored by",strsplit(formula," ")[[1]][1]))

  form<-as.formula(paste("t(c) ~", formula))
  mod<-vegan::cca(form, data=i)
  plot(mod, dis="bp", scaling="sites", type="n")
  points(mod, pch=21, col=rgb(0,0,0,.05), cex=.4, "sp",...)
  text(mod, dis="bp")
  return(mod)
}
