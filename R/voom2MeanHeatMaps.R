#' @title PCA of voom-transformed counts.
#'
#'
#' @description
#' \code{voom2PCA} Run principal component analysis on matrix of voom-normalized counts
#'
#' @param v Counts matrix, typically transformed by limma::voom.
#' Possibly output from pipelimma, in the slot voom[["E"]].
#' @param grps The grouping factor within which to calculate mean expression values
#' @param thresh The threshold at which to cut the dendrogram into
#' separate groups. This also defines the number of sub-heatmaps to make.
#' Smaller values lead to more groups and more sub heatmaps.
#' @param calcMeans Logical, should means for each group be calculated?
#' Simplifies plotting considerably.
#' @param newIDs If calcMeans=FALSE, a vector of names for the columns in the heatmaps
#' @param allHMtoo Should the complete heatmap be plotted in addition to the sub heatmaps?
#'
#' @details This function is typically used to calculate group mean
#' normalized expression to display the effects of different
#' experimental factors. Specifically, for each level it calculates the mean.
#' Then normalizes across levels using scale(..., center=TRUE, scale=TRUE).
#' These scaled, centered values are used to calculate a dendrogram with hclust.
#' The dendrogram is split into groups using cuttree. The dendrogram is plotted,
#' colored by groups. Complete and sub heatmaps are subsequently plotted.

#' @return a dataframe with the normalized means.
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' set.seed(42)
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' ### Not Run ### stats.fullmodel<-stats$stats
#' ### Not Run ### which.toplot<-which(stats.fullmodel$ebayesQvalue_treatmentTumor<=1e-20)
#' ### Not Run ### v<-stats$voom[["E"]]
#' ### Not Run ### v.means<-voom2MeanHeatMaps(v=v[which.toplot,], grps=info$rep, thresh=21)
#' @importFrom  gdata read.xls startsWith
#' @importFrom  Heatplus annHeatmap2
#' @export
voom2MeanHeatMaps<-function(v, grps=info$Treatment, calcMeans=T){
  opar<-par(no.readonly=T)
  if(calcMeans){
    v<-as.matrix(v)
    colnames(v)[startsWith(colnames(v),"X")]<-gsub("X","",colnames(v)[startsWith(colnames(v),"X")])
    mean.v<-v[,0]
    for(i in unique(grps)){
      temp<-data.frame(rowMeans(v[,grps==i]))
      colnames(temp)<-i
      mean.v<-data.frame(mean.v,temp)
    }
  }else{
    mean.v<-v
  }
  scl<-t(mean.v)
  scl2<-scale(scl, center=T, scale=F)
  std.mean.v<-t(scl2)
  pal<-colorRampPalette(c('dark blue','white','dark red'))
  v2<-data.matrix(std.mean.v)
  rownames(v2)<-NULL
  hm1<-annHeatmap2(x=v2, labels=list(row=F), scale="none",col=pal, legend=T)
  plot(hm1)
  par(opar)
}
