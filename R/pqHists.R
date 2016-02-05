#' @title P- and q-value distributions plotted side-by-side
#'
#'
#' @description
#' \code{pqHists} Make histograms of p-value and q-value distributions
#'
#' @param x matrix containing p-values and q-values
#' @param what.p the name of the column containing p-values
#' @param what.q the name of the column containing q-values
#' @param main The stem of the title, to which significance counts are attached
#' @param ... additional arguments passed on to hist, like the color of the bars
#'
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' ### Not Run ### stats.fullmodel<-stats$stats
#' ### Not Run ### pqHists(stats.fullmodel, what.p="ebayesPvalue_treatmentTumor",
#' ### Not Run ### what.q="ebayesQvalue_treatmentTumor", main="main effect treatment", breaks=100)
#'
#' @export
pqHists<-function(x, what.p="p.value",what.q="q.value",alpha=0.05,main="significance distribution",...){
  par(mfrow=c(1,2))
  hist(x[,grep(what.p, colnames(x))],
       xlab="p.value",
       main=paste(main,
                  "\n","n.pvalue < ",
                  alpha," = ",
                  sum(x[,grep(what.p, colnames(x))]<=alpha),sep=""),
       ...)
  hist(x[,grep(what.q, colnames(x))],
       xlab="q.value",
       main=paste(main,
                  "\n","n.qvalue < ",
                  alpha," = ",
                  sum(x[,grep(what.q, colnames(x))]<=alpha),sep=""),
       ...)
  par(mfrow=c(1,1))
}
