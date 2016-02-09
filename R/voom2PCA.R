#' @title PCA of voom-transformed counts.
#'
#'
#' @description
#' \code{voom2PCA} Run principal component analysis on matrix of voom-normalized counts
#'
#' @param v Counts matrix, typically transformed by limma::voom.
#' Possibly output from pipelimma, in the slot "voom".
#' @param info The experimental design information matrix
#' @param ids A vector of the individual names
#' @param plotit Logical, should the pca be plotted?
#' @param ... additional arguments passed on to pairs, for example, the colors of bars
#'
#' @details This function uses the R function princomp to calculate principal components

#' @return a list containing a dataframe with the experimental design data, merged with the 1st
#' 3 principal component axes and a vector of %variance explained by PCA axes.
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' ### Not Run ### pc <- voom2PCA(v=stats$voom[["E"]], info=info, ids=rownames(info),plotit=TRUE)
#' ### Not Run ### library(ggplot2)
#' ### Not Run ### ggplot(pc, aes(x=PC1, y=PC2, col=treatment))+geom_point()
#' @export
voom2PCA<-function(v, info, ids, plotit=TRUE,pcas2return=3,plot.cols="black",...){

  pc<-prcomp(t(v))
  prop.var<-round(((pc$sdev)^2 / sum(pc$sdev^2))*100,1)
  dat<-data.frame(info, pc$x[,1:pcas2return])
  if(plotit){
    par(mfrow=c(2,2))
#     bp<-barplot(prop.var[1:5], ylab="% Variance Explained",
#                 main="distribution of PCA effects", xlab="PCA Axis")
#     axis(1, at=bp[,1], labels=1:5, title)
#     with(dat, plot(PC1,PC2, type="n", bty="n", main="PC1 vs. PC2"))
#     with(dat, text(PC1,PC2, label=ids))
#     with(dat, plot(PC1,PC3, type="n", bty="n", main="PC1 vs. PC3"))
#     with(dat, text(PC1,PC3, label=ids))
#     with(dat, plot(PC2,PC3, type="n", bty="n", main="PC2 vs. PC3"))
#     with(dat, text(PC2,PC3, label=ids))
    par(mfrow=c(1,1))
    pairs(dat[,grepl("PC", colnames(dat))], col=plot.cols,
          labels=paste(colnames(dat)[grepl("PC", colnames(dat))]," ",prop.var,"%",sep=""),...)
  }

  return(list(dat,prop.var))
}
