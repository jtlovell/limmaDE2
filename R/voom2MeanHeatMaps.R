#' @title PCA of voom-transformed counts.
#'
#'
#' @description
#' \code{voom2PCA} Run principal component analysis on matrix of voom-normalized counts
#'
#' @param v Counts matrix, typically transformed by limma::voom. Possibly output from pipelimma, in the slot voom[["E"]].
#' @param grps The grouping factor within which to calculate mean expression values
#' @param thresh The threshold at which to cut the dendrogram into separate groups. This also defines the number of sub-heatmaps to make.
#' Smaller values lead to more groups and more sub heatmaps.
#' @param calcMeans Logical, should means for each group be calculated? Simplifies plotting considerably.
#' @param newIDs If calcMeans=FALSE, a vector of names for the columns in the heatmaps
#' @param allHMtoo Should the complete heatmap be plotted in addition to the sub heatmaps?
#'
#' @details This function is typically used to calculate group mean normalized expression to display the effects of different
#' experimental factors. Specifically, for each level it calculates the mean. Then normalizes across levels using
#' scale(..., center=TRUE, scale=TRUE). These scaled, centered values are used to calculate a dendrogram with hclust.
#' The dendrogram is split into groups using cuttree. The dendrogram is plotted, colored by groups. Complete and sub heatmaps
#' are subsequently plotted.

#' @return a dataframe with the normalized means.
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' stats.fullmodel<-stats$stats
#' which.toplot<-which(stats.fullmodel$ebayes_treatmenttumor_q.value<=1e-20)
#' v<-stats$voom[["E"]]
#' v.means<-voom2MeanHeatMaps(v=v[which.toplot,], grps=info$rep, thresh=21)
#' @export
voom2MeanHeatMaps<-function(v, grps=info$Treatment, thresh=7, calcMeans=T, newIDs=NA, allHMtoo=T){
  opar<-par()

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
    colnames(mean.v)<-newIDs
    mean.v<-mean.v[,order(newIDs)]
  }
  scl<-t(mean.v)
  scl2<-scale(scl, center=T, scale=F)
  std.mean.v<-t(scl2)
  hc<-hclust(dist(std.mean.v))
  d.out=colour_clusters(hc,h=thresh)
  cut.in<-cutree(hc,h=thresh)
  par(mfrow=c(1,1))
  pal<-colorRampPalette(c('dark blue','white','dark red'))
  plot(d.out)
  abline(h=thresh, lty=2)
  if(allHMtoo){
    v2<-data.matrix(std.mean.v)
    rownames(v2)<-NULL
    hm1<-annHeatmap2(x=v2, labels=list(row=F), scale="none",col=pal, legend=T)
    plot(hm1, cex=.5)
  }
  for(i in names(table(cut.in))[table(cut.in)>1]){
    temp<-std.mean.v[cut.in==i,]
    rownames(temp)<-NULL
    hm1<-annHeatmap2(x=data.matrix(temp), labels=list(row=F), scale="none",col=pal, legend=T)
    plot(hm1, cex.axis=.5, main = paste("cluster",i), col.main=rainbow(length(names(table(cut.in))[table(cut.in)>1]))[i])
  }
  return(mean.v)
  par(opar)
}
