#' @title A pipeline for LIMMA.
#'
#'
#' @description
#' \code{pipeLIMMA} Run a pipeline of LIMMA functions for differential gene expression.
#'
#' @param counts A count matrix
#' @param info An experimental design matrix
#' @param formula A character string that can be coerced to a formula. Specify if a contrast
#' model is not desired.
#' @param design A design matrix, usually created by a call from "model.matrix". Used only if
#' limma:contrast.fit is intended as the statistical modeling function. NULL values force
#' a traditional test of effects via limma::lmFit/ebayes.
#' @param contrast.matrix A matrix of contrasts, usually created by a call
#' from limma:makeContrasts. Used only limma:contrast.fit is intended as the statistical
#' modeling function. NULL values force a traditional test of effects via limma::lmFit/ebayes.
#' @param block A string that represents an individual that was repeatedly measured,
#' if NULL, runs the analysis without a blocking / duplicate correlation factor
#' @param use.qualityWeights Logical, run voom with quality weights or not?
#' @param geneIDs The names of genes. If NA, use row names from counts matrix
#' @param getTopTable Logical, return toptable statistics?
#' @param getEbayes Logical, return ebayes statistics?
#' @param simplify Logical, return a element with the F-statistics from the main model?
#' @param verbose Logical, return progress updates?
#' @param ... additional arguments, not currently in use.
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{calculate normalization factors via edgeR::calcNormFactors}
##'  \item{2. }{Run limma::voom transformation}
##'  \item{3. }{Run limma:lmFit linear modeling}
##'  \item{4. }{Run limma::ebayes statistical modeling}
##'  \item{5. }{Collect Log2 Fold-Changes using limma:topTable}
##'  \item{6. }{Collate and ouput statistics and voom transformed data}
##' }
##'
#' @return a list with 2 elements (if simple=TRUE)
##' \itemize{
##'  \item{"stats"}{: the statsistics generated from ebayes and topTable}
##'  \item{"voom"}{: the voom normalized counts data}
##' }
##'
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=info$rep)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' @importFrom  edgeR calcNormFactors DGEList
#' @importFrom  qvalue qvalue
#' @export
pipeLIMMA<-function(counts, info, formula=NULL, contrast.matrix=NULL, block=NULL,
                    design=NULL, use.qualityWeights=TRUE,
                    geneIDs=NA, getTopTable=FALSE, getEbayes=TRUE,
                    simplify=TRUE, verbose=TRUE, ...){

  if(is.null(block)) {
    useBlock=FALSE
  }else{
    useBlock=TRUE
  }

  if(verbose) cat("calculating normalization factors ... \n")
  if(is.na(geneIDs)){
    geneIDs<-rownames(counts)
  }
  if(is.null(design)){
    design<-model.matrix(as.formula(formula), data = info)
  }
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)

  if(use.qualityWeights){
    if(verbose) cat("running voom normalization correcting for quality weights ... \n")
    v <- voomWithQualityWeights(y, design=design, plot = T)
  }else{
    if(verbose) cat("running voom normalization ... \n")
    v <- voom(y, design=design, plot = T)
  }
  if(useBlock){
    if(verbose) cat("calculating duplicate correlation among replicates ... \n")
    dupcor <- duplicateCorrelation(counts,design, block=as.factor(block))
    if(!is.null(contrast.matrix)){
      if(verbose) cat("fitting model to contrast matrix ... \n")
      fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block))
      fit <- contrasts.fit(fit, contrast.matrix)
      fit <- eBayes(fit)
    }else{
      if(verbose) cat("fitting linear model ... \n")
      fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block))
      fit <- eBayes(fit[,-1])
    }
  }else{
    if(!is.null(contrast.matrix)){
      if(verbose) cat("fitting model to contrast matrix ... \n")
      fit <- lmFit(v, design=design)
      fit <- contrasts.fit(fit, contrast.matrix)
      fit <- eBayes(fit)
    }else{
      if(verbose) cat("fitting linear model ... \n")
      fit <- lmFit(v, design=design)
      fit <- eBayes(fit[,-1])
    }
  }
  if(verbose) cat("processing statistics and calculating q-values ... \n")

  main.out<-data.frame(gene=geneIDs,
                  sigma=fit$sigma,
                  s2.post=fit$s2.post,
                  Amean=fit$Amean,
                  Fstat=fit$F,
                  Fpvalue=fit$F.p.value,
                  Fqvalue=p.adjust(fit$F.p.value, method="BH"))

  ebayes.coef<-fit$coefficients
  colnames(ebayes.coef)<-paste("ebayesCoef_",colnames(ebayes.coef),sep="")

  ebayes.t<-fit$t
  colnames(ebayes.t)<-paste("ebayesTstat_",colnames(ebayes.t),sep="")

  ebayes.p<-fit$p.value
  colnames(ebayes.p)<-paste("ebayesPvalue_",colnames(ebayes.p),sep="")

  ebayes.q<-apply(ebayes.p, 2, function(x) p.adjust(x, method="BH"))
  colnames(ebayes.q)<-gsub("ebayesPvalue_","ebayesQvalue_",colnames(ebayes.p))

  coefnames<-colnames(fit)
  lfcs<-lapply(coefnames, function(x) {
    tt<-topTable(fit, coef=x, p.value=1, number=100000)
    tt<-tt[,c("logFC","AveExpr")]
    colnames(tt)<-paste(x,colnames(tt),sep="_")
    tt<-tt[match(geneIDs, row.names(tt)),]
  })
  lfcs<-do.call(cbind, lfcs)

  all.out<-data.frame(main.out, lfcs, ebayes.coef, ebayes.t, ebayes.p, ebayes.q)

  return(list(stats=all.out, voom=v))
}
