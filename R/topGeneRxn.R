#' @title Plots of expression vs. phenotype data.
#'
#' @description
#' \code{volcanoPair} Make a matrix of plots for a set of genes. Typically,
#' for a given effect, pull out the strongest genes and plot their phenotype vs.
#' voom normalized gene expression data. Current implementation works best for
#' continuous x data.
#'
#' @param v Voom-transformed (or similar) expression data
#' @param info Experimental design dataframe or matrix
#' @param sig A vector of transformed (or not) p-values used to define
#' significance
#' @param nsig The number of genes to plot
#' @param alpha A threshold for significance, currently not used.
#' @param xdat The column in the info dataset to plot on the x axis
#' @param coldat The column in the info dataset to color by, should be a factor
#' @param pointcols the colors for the factor (coldat). Needs to match exactly the number
#' of levels in coldat
#' @param ylab The label for the y-axis
#' @param xlab The label for the x-axis
#' @param title The title for the plot
#' @param ... additional arguments passed on to geom_point()
#'
#' @details Takes a gene expression matrix, subsets by the nsig most significant genes,
#' melts to the long format, then plots via ggplot in a scatterplot matrix.

#' @return a table with the data used to make the plot
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment + rep", block=NULL, getTopTable=T, getEbayes=T)
#' stats<-stats$stats
#'

#' @import  ggplot
#' @importFrom  reshape2 melt
#' @export
topGeneRxn<-function(v, info, sig, xdat, coldat, alpha=0.05, nsig=20, geneIDs=NULL,
                     pointcols=c("darkred","forestgreen","cornflowerblue"),
                     ylab="voom normalized expression",
                     xlab="phenotype", title="comparison of expression and phenotype",...){

  vs<-v
  if(!is.null(geneIDs)){
    vs<-vs[geneIDs,]
  }else{
    best<-order(sig)[1:nsig]
    vs<-vs[best,]
  }

  df<-cbind(info[,c(xdat,coldat)], t(vs))
  vtp<-melt(df, id.vars=c(xdat,coldat))
  print(
    ggplot(vtp, aes_string(x=xdat, y="value"))+ geom_point(aes_string(col=coldat))+
      facet_wrap(~variable, scales="free_y", nrow=5, ncol=4)+
      scale_color_manual(values=pointcols)+
      stat_smooth(span = 200,se=F, lty=2, alpha=.2, lwd=.5, col="black")+
      theme_bw()+
      theme()+
      scale_y_continuous(ylab)+
      scale_x_continuous(xlab)+
      theme(panel.grid.major = element_blank() ,
            panel.grid.minor = element_blank(),
            strip.text.x = element_text(size = 8),
            strip.background = element_blank())+
      ggtitle(title)
  )
  return(vtp)
}
