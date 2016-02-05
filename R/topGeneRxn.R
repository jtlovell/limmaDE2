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
#' @param paletteChoice If using a palette instead of predefined colors, which r color
#' brewer palette to use?
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
#' @import  ggplot
#' @importFrom  reshape2 melt
#' @export
topGeneRxn<-function(v, info, sig, xdat=NULL, coldat=NULL, alpha=0.05, nsig=20, geneIDs=NULL,paletteChoice=NULL,
                     pointcols=c("darkred","forestgreen","cornflowerblue"),
                     ylab="voom normalized expression",
                     xlab=NULL, title="comparison of expression and phenotype",...){
  if(is.null(xlab)) xlab<-xdat
  if(is.null(xdat) & is.null(coldat)) cat("either xdat or coldat must be provided \n")
  if(!is.null(pointcols) & !is.null(paletteChoice)) cat("both manual and palette colors requests, using palette \n")

  if(is.null(coldat)){
    ncols<-length(unique(info[,xdat]))
  }else{
    ncols<-length(unique(info[,coldat]))
  }

  if(!is.null(paletteChoice)){
    pointcols<-brewer.pal(n=ncols, name=paletteChoice)
  }

  if(sum(is.na(info[,xdat]))>0 | sum(is.na(info[,coldat]))>0){
    ind<-which(complete.cases(info[,c(xdat,coldat)]))
    vs<-v[,ind]
    info<-info[ind,]
  }else{
    vs<-v
  }

  if(!is.null(geneIDs)){
    vs<-vs[geneIDs,]
  }else{
    best<-order(sig)[1:nsig]
    vs<-vs[best,]
  }

  if(identical(info[,coldat],info[,xdat]) | is.null(xdat) | is.null(coldat)){
    df<-data.frame(as.character(info[,c(xdat)]), t(vs))
    colnames(df)[1]<-xdat
    vtp<-melt(df, id.vars=c(xdat))
    vtp$value<-as.numeric(as.character(vtp$value))
    print(
      ggplot(vtp, aes_string(x=xdat, y="value"))+ geom_boxplot(aes_string(col=xdat))+
        facet_wrap(~variable, scales="free_y", nrow=5, ncol=4)+
        scale_color_manual(values=pointcols)+
        stat_smooth(span = 200,se=F, lty=2, alpha=.2, lwd=.5, col="black")+
        theme_bw()+
        labs(x = xlab, y=ylab, title=title)+
        theme(panel.grid.major = element_blank() ,
              panel.grid.minor = element_blank(),
              strip.text.x = element_text(size = 8),
              strip.background = element_blank(),
              axis.text.x = element_text(angle = 90,vjust=.5, hjust=1))
    )
  }else{
    df<-data.frame(info[,c(xdat, coldat)], t(vs))
    vtp<-melt(df, id.vars=c(xdat,coldat))
    vtp$value<-as.numeric(as.character(vtp$value))
    print(
      ggplot(vtp, aes_string(x=xdat, y="value"))+ geom_point(aes_string(col=coldat))+
        facet_wrap(~variable, scales="free_y", nrow=5, ncol=4)+
        scale_color_manual(values=pointcols)+
        stat_smooth(span = 200,se=F, lty=2, alpha=.2, lwd=.5, col="black")+
        theme_bw()+
        labs(x = xlab, y=ylab, title=title)+
        theme(panel.grid.major = element_blank() ,
              panel.grid.minor = element_blank(),
              strip.text.x = element_text(size = 8),
              strip.background = element_blank())
    )
  }
  return(vtp)
}
