#' @title Volcano plot.
#'
#' @description
#' \code{volcanoPlot} Make an expression volcano plot
#' from log2-fold change and p-value data
#'
#' @param pval A vector of p-values
#' @param lfc A vector of log-2 fold change (or similar data)
#' @param sig Either a 0/1 binary vector of significance, or a
#' vector of transformed p-values
#' @param alpha If providing a vector of transformed p-values,
#' this specifies the threshold for siginficance. Otherwise, not used.
#' @param pointcols Vector of legnth 2. The first element indicates
#' the color of significant points and the second the non sigificant points
#' @param pchs Follows pointcols for point shapes
#' @param cex Follows pointcols for point sizes
#' @param legpos The position of the legend. Defaults to topright
#' @param leginset How much the legend should be inset - defaults to just
#' outside the plotting window.
#' @param ylab The label for the y axis
#' @param xlab The label for the x axis
#' @param ... additional arguments passed on to plot
#'
#' @details Using base R graphics, this function tabulates
#' the number of significant values and plots

#' @return a table with the number and direction of significance
#' @examples
#' \dontrun{
#' data(kidney) #from the simseq package
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info,
#'                  formula = " ~ treatment",
#'                  block=NULL,
#'                  getTopTable=T,
#'                  getEbayes=T)
#' stats<-stats$stats
#' volcanoPlot(pval=stats$ebayesPvalue_treatmentTumor,
#'             lfc=stats$treatmentTumor_logFC,
#'             sig=stats$ebayesQvalue_treatmentTumor,
#'             alpha=0.0005, pointcols=c("blue","grey"), pchs=c(4,6))
#' }
#' @export
volcanoPlot<-function(pval, lfc, sig, alpha=0.05,
                      pointcols=c(rgb(1,0,0,.5),rgb(0,0,0,.5)),
                      pchs=c(19,19), cex=c(.5,.5),
                      legpos="topright",leginset = c(-0.25,0),
                      ylab="-log10 P-value", xlab="log2 Fold Change",...){

  if(length(unique(sig))>2){
    sig<-ifelse(sig<=alpha, TRUE, FALSE)
  }else{
    sig<-ifelse(sig==1, TRUE, FALSE)
  }
  updown<-ifelse(lfc<=0,"down","up")
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  udsig<-table(updown[sig])
  ns<-sum(!sig)
  plot(x=lfc, y=-log10(pval), type="n", ylab=ylab, xlab=xlab, ...)
  points(x=lfc[sig], y=(-log10(pval))[sig], col=pointcols[1], pch=pchs[1], cex=cex[1])
  points(x=lfc[!sig], y=(-log10(pval))[!sig], col=pointcols[2], pch=pchs[2], cex=cex[2])
  legend(legpos, inset=leginset[1:2], cex=.8,
         legend=c(paste("< 0 : ", udsig[1], sep=""),
                  paste("> 0 : ", udsig[2], sep=""),
                  paste("NS : ", ns, sep="")),
         pch=c(pchs[1],pchs[1],pchs[2]),
         col=c(pointcols[1], pointcols[1], pointcols[2]),
         ncol=1, title="significance", xjust=0, yjust=.5, bty="n")

  return(table(updown,sig))
}
