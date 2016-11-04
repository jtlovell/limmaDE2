#' @title Scatter plot of log2 fold changes.
#'
#' @description
#' \code{volcanoPair} Make an expression volcano plot
#' from log2-fold change and p-value data
#'
#' @param lfc1 A vector of log-2 fold change (or similar data) for the x axis
#' @param lfc2 A vector of log-2 fold change (or similar data) for the y axis
#' @param cols A vector of colors,
#' if specified, overrides all other coloring parameters
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
#' @details Using base R graphics, this function tabulates
#' the number of significant values and plots

#' @return a table with the number of significant effects
#' @examples
#' \dontrun{
#' data(kidney) #from the simseq package
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info,
#'                  formula = " ~ treatment + rep",
#'                  block=NULL, getTopTable=T, getEbayes=T)
#' stats<-stats$stats
#' volcanoPair(lfc1=stats$treatmentTumor_logFC,
#'            lfc2=stats$rep6090_logFC,
#'             sig1=stats$ebayesQvalue_treatmentTumor,
#'             sig2=stats$ebayesQvalue_rep6090)
#' }
#'
#' @export
volcanoPair<-function(lfc1, lfc2,
                      legpos="topleft",leginset=c(0,0), legcex=1,
                      pt.col = NULL, pt.pch = NULL, pt.cex = NULL,
                      line.col = "grey", line.lty =2, line.lwd = 1, ...){
  plot(x=lfc1, y=lfc2, type="n",...)
  lines(c(0,0),c(min(lfc2),max(lfc2)), lty=line.lty, col=line.col, lwd = line.lwd)
  lines(c(min(lfc1),max(lfc1)),c(0,0), lty=line.lty, col=line.col, lwd = line.lwd)

  if(is.null(pt.col)) pt.col<-"grey"
  if(is.null(pt.pch)) pt.pch<-19
  if(is.null(pt.cex)) pt.cex<-1

  points(x=lfc1, y=lfc2, col = pt.col, pch = pt.pch, cex = pt.cex)

  if(sum(c(length(pt.col),length(pt.pch),length(pt.cex)))>1){
    if(any(all(length(pt.col) > 1, length(pt.pch) > 1),
           all(length(pt.col) > 1, length(pt.cex) > 1))){
      tab<-table(pt.col)
    }else{
      if(length(pt.col)>1) {
        tab<-table(pt.col)
        leg.col<-names(tab)
        leg.pch = pt.pch
        leg.cex = pt.cex
      }else{
        if(length(pt.pch)>1){
          tab<-table(pt.pch)
          leg.pch<-names(tab)
          leg.col = pt.col
          leg.cex = pt.cex
        }else{
          if(length(pt.cex)>1){
            tab<-table(pt.cex)
            leg.cex<-names(tab)
            leg.col = pt.col
            leg.pch = pt.pch
          }
        }
      }
    }

    legend(legpos, leginset, tab,
           col = leg.col,
           pch = leg.pch,
           pt.cex = leg.cex,
           bty = "n")
    return(tab)
  }
}
