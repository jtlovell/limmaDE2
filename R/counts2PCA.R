#' @title PCA of voom-transformed counts.
#'
#' @description
#' \code{counts2PCA} Run principal component analysis on matrix of
#' voom-normalized counts
#'
#' @param counts Counts matrix, typically transformed by limma::voom.
#' Possibly output from pipelimma, in the slot "voom".
#' @param info The experimental design information matrix
#' @param ids A vector of the individual names
#' @param plotit Logical, should the pca be plotted?
#' @param ... additional arguments passed on to pairs, for example,
#' the colors of bars
#'
#' @details This function uses the R function princomp to calculate
#' principal components

#' @return a list containing a dataframe with the experimental
#' design data, merged with the 1st 3 principal component axes
#' and a vector of %variance explained by PCA axes.
#' @examples
#' \dontrun{
#' data(kidney) #from the simseq package
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic,
#'                  treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts,
#'                  info=info,
#'                  formula = " ~ treatment",
#'                  block=NULL)
#' pc <- voom2PCA(counts=stats$voom[["E"]],
#'                info=info,
#'                ids=rownames(info),
#'                plotit=TRUE)
#' library(ggplot2)
#' ggplot(pc, aes(x=PC1, y=PC2, col=treatment))+
#'    geom_point()
#' }
#' @export
counts2PCA<-function(counts, info, ids, plotit=TRUE,
                     pcas2return=3,plot.cols="black",...){

  pc<-prcomp(t(counts))
  prop.var<-round(((pc$sdev)^2 / sum(pc$sdev^2))*100,1)
  dat<-data.frame(info, pc$x[,1:pcas2return])
  if(plotit){
    pairs(dat[,grepl("PC", colnames(dat))], col=plot.cols,
          labels=paste(colnames(dat)[grepl("PC", colnames(dat))],
                       " ",
                       prop.var,"%",sep=""),...)
  }
  return(list(dat,prop.var))
}
