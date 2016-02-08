#' @title Scatter plot of log2 fold changes.
#'
#' @description
#' \code{volcanoPair} Make an expression volcano plot from log2-fold change and p-value data
#'
#' @param lfc1 A vector of log-2 fold change (or similar data) for the x axis
#' @param lfc2 A vector of log-2 fold change (or similar data) for the y axis
#' @param sig1 Either a 0/1 binary vector of significance, or a
#' vector of transformed p-values for coloring the x axis points
#' @param sig2 Either a 0/1 binary vector of significance, or a
#' vector of transformed p-values for coloring the y axis points
#' @param alpha If providing a vector of transformed p-values,
#' this specifies the threshold for siginficance. Otherwise, not used.
#' @param pointcols Vector of legnth 4. Indicates  the color of significant points
#' for both [1], just the y axis [2], just the x axis [3] or neighter.
#' @param pchs Follows pointcols for point shapes
#' @param cex Follows pointcols for point sizes
#' @param legpos The position of the legend. Defaults to the top right.
#' @param ... additional arguments passed on to plot
#'
#' @details Using base R graphics, this function tabulates the number of significant values and plots

#' @return a table with the number of significant effects
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment + rep", block=NULL, getTopTable=T, getEbayes=T)
#' ### Not Run ### stats<-stats$stats
#' ### Not Run ### volcanoPair(lfc1=stats$treatmentTumor_logFC,
#' ### Not Run ###             lfc2=stats$rep6090_logFC,
#' ### Not Run ###             sig1=stats$ebayesQvalue_treatmentTumor,
#' ### Not Run ###             sig2=stats$ebayesQvalue_rep6090)
#' @export
volcanoPair<-function(lfc1, lfc2, sig1, sig2, alpha=0.05,
                      pointcols=c("darkred","skyblue","forestgreen",rgb(0,0,0,.5)),
                      pchs=c(19,19,19,19), cex=c(.5,.5,.5,.5),legpos="topleft"
                      ,...){
  if(length(unique(c(sig1,sig2)))>2){
    sig1<-ifelse(sig1<=alpha, TRUE, FALSE)
    sig2<-ifelse(sig2<=alpha, TRUE, FALSE)
  }else{
    sig1<-ifelse(sig1==1, TRUE, FALSE)
    sig2<-ifelse(sig2==1, TRUE, FALSE)
  }
  sigs<-ifelse(sig1 & sig2, "both",
               ifelse(sig1,"sig1",
                      ifelse(sig2,"sig2","NS")))
  tab<-table(factor(sigs,levels=c("both","sig1","sig2","NS")))
  plot(x=lfc1, y=lfc2, type="n", ...)
  lines(c(0,0),c(min(lfc2),max(lfc2)), lty=2, col="grey")
  lines(c(min(lfc1),max(lfc1)),c(0,0), lty=2, col="grey")
  points(x=lfc1[sigs=="NS"], y=lfc2[sigs=="NS"], col=pointcols[4], pch=pchs[4], cex=cex[4])
  points(x=lfc1[sigs=="sig2"], y=lfc2[sigs=="sig2"], col=pointcols[3], pch=pchs[3], cex=cex[3])
  points(x=lfc1[sigs=="sig1"], y=lfc2[sigs=="sig1"], col=pointcols[2], pch=pchs[2], cex=cex[2])
  points(x=lfc1[sigs=="both"], y=lfc2[sigs=="both"], col=pointcols[1], pch=pchs[1], cex=cex[1])

  legend(legpos, cex=.8,
         legend=c(tab[1], tab[2], tab[3], tab[4]),
         pch=c(pchs[1],pchs[2],pchs[3], pchs[4]),
         col=c(pointcols[1], pointcols[2], pointcols[3], pointcols[4]),
         ncol=2, title="n genes", xjust=0, yjust=.5, bty="n")

  return(tab)
}
