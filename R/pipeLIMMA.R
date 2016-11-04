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
#' from limma:makeContrasts. Used only when limma:contrast.fit is intended as the statistical
#' modeling function. NULL values force a traditional test of effects via limma::lmFit/ebayes.
#' @param block A string that represents an individual that was repeatedly measured,
#' if NULL, runs the analysis without a blocking / duplicate correlation factor
#' @param runVoom Logical, if TRUE, normalizes the counts matrix via voom. If FALSE,
#' assumes the counts matrix is already voom-normalized. Pre-running voom will speed up
#' analyses with multiple pipeLIMMA calls.
#' @param use.topTable Logical, report F-statistics across all factors? If true,
#' a third element is returned called fstats.
#' @param use.qualityWeights Logical, run voom with quality weights or not?
#' @param geneIDs The names of genes. If NA, use row names from counts matrix
#' @param verbose Logical, return progress updates?
#' @param plotVoom Logical, plot the voom fit? Defaults to FALSE
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
##'  \item{"fstats"}{: if a formula is provided,
##'  the toptable returned fstatsistics across all estimates of each factor}
##' }
##'
#' @examples
#' \dontrun{
#' data(kidney) # from simseq
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic,
#'                  treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts,
#'                  info=info,
#'                  formula = " ~ treatment",
#'                  block=info$rep)
#' stats<-pipeLIMMA(counts=counts,
#'                  info=info,
#'                  formula = " ~ treatment",
#'                  block=NULL)
#' }
#' @importFrom  edgeR calcNormFactors DGEList
#' @importFrom  qvalue qvalue
#' @export
pipeLIMMA<-function(counts, info, formula=NULL,
                    contrast.matrix=NULL, block=NULL,
                    design=NULL, runVoom=TRUE,
                    use.qualityWeights=TRUE,use.topTable=FALSE,
                    geneIDs=NA, verbose=TRUE, plotVoom=FALSE, ...){

  extractTopTable<-function(fit, formula){
    n.effects<-nchar(formula)-nchar(gsub("+","",formula, fixed=T))+1
    effectNames<-sapply(1:n.effects, function(x) {
      out<-strsplit(formula,"[+]")[[1]][x]
      out<-gsub(" ","",out, fixed=T)
      out<-gsub("~","",out, fixed=T)
      out<-gsub("*","_x_",out, fixed=T)
      out})
    maxTermInteraction<-max(sapply(effectNames, function(x)
      (nchar(x) - nchar(gsub("_x_","",x)))/3))

    if(maxTermInteraction>1)
      stop("use.topTable is not possible with greater than a 2-way interaction")
    colids<-colnames(fit$p.value)
    colidInt<-colids[grep(":",colids)]
    colidMain<-colids[grep(":",colids, invert=T)]
    ttColNames<-c("logFC",  "AveExpr",  "t",  "P.Value",  "adj.P.Val", "B", "F")
    tt<-lapply(effectNames, function(i){
      if(grepl("_x_", i)){
        temp<-sapply(1:2, function(x) strsplit(i,"_x_")[[1]][x])
        toget<-colidInt[grepl(temp[1], colidInt) & grepl(temp[2], colidInt)]
        wh<-which(colids %in% toget)
        tt<-topTableF(fit[,wh],sort="none",n=Inf,adjust.method="none")
      }else{
        wh<-grep(i, colidMain)
        tt<-topTableF(fit[,wh],sort="none",n=Inf,adjust.method="none")
      }
      tt<-tt[,colnames(tt) %in% ttColNames]
      tt$adj.P.Val<-NULL
      tt$Q.Value<-p.adjust(tt$P.Value, method = "fdr")
      colnames(tt)<-paste(i, colnames(tt),sep="_")
      tt$gene<-row.names(tt)
      return(tt)
    })
    if(length(tt)>1){
      out<-merge(tt[[1]],tt[[2]],by="gene")
      if(length(tt)>2){
        for(i in 3:length(tt)){
          out<-merge(out,tt[[i]],by="gene")
        }
      }
    }else{
      out<-tt[[1]]
    }
    return(out)
  }

  if(is.null(block)) {
    useBlock=FALSE
  }else{
    useBlock=TRUE
  }
  if(is.na(geneIDs)){
    geneIDs<-rownames(counts)
  }
  if(is.null(design)){
    design<-model.matrix(as.formula(formula), data = info)
  }
  if(runVoom){
    if(verbose) cat("calculating normalization factors ... \n")
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y)
    if(use.qualityWeights){
      if(verbose) cat("running voom normalization correcting for quality weights ... \n")
      v <- voomWithQualityWeights(y, design=design, plot = plotVoom)
    }else{
      if(verbose) cat("running voom normalization ... \n")
      v <- voom(y, design=design, plot = plotVoom)
    }
  }else{
    v<-counts
  }

  if(useBlock){
    if(verbose) cat("calculating duplicate correlation among replicates ... \n")
    dupcor <- duplicateCorrelation(counts,design, block=as.factor(block))
    if(!is.null(contrast.matrix)){
      if(verbose) cat("fitting model to contrast matrix ... \n")
      fit <- lmFit(v,
                   design=design,
                   correlation=dupcor$consensus,
                   block=as.factor(block), ...)
      fit <- contrasts.fit(fit, contrast.matrix)
      fit <- eBayes(fit)
    }else{
      if(verbose) cat("fitting linear model ... \n")
      fit <- lmFit(v,
                   design=design,
                   correlation=dupcor$consensus,
                   block=as.factor(block), ...)
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
                       Fqvalue=p.adjust(fit$F.p.value, method="fdr"))

  ebayes.coef<-fit$coefficients
  colnames(ebayes.coef)<-paste("ebayesCoef_",colnames(ebayes.coef),sep="")

  ebayes.t<-fit$t
  colnames(ebayes.t)<-paste("ebayesTstat_",colnames(ebayes.t),sep="")

  ebayes.p<-fit$p.value
  colnames(ebayes.p)<-paste("ebayesPvalue_",colnames(ebayes.p),sep="")

  ebayes.q<-apply(ebayes.p, 2, function(x) p.adjust(x, method = "fdr"))
  colnames(ebayes.q)<-gsub("ebayesPvalue_","ebayesQvalue_",colnames(ebayes.p))

  coefnames<-colnames(fit)
  lfcs<-lapply(coefnames, function(x) {
    tt<-topTable(fit, coef=x,sort="none",n=Inf)
    tt<-tt[,c("logFC","AveExpr")]
    colnames(tt)<-paste(x,colnames(tt),sep="_")
    tt<-tt[match(geneIDs, row.names(tt)),]
  })
  lfcs<-do.call(cbind, lfcs)

  all.out<-data.frame(main.out, lfcs, ebayes.coef, ebayes.t, ebayes.p, ebayes.q)

  if(is.null(formula)){
    return(list(stats=all.out, voom=v))
  }else{
    fstats<-extractTopTable(fit=fit, formula=formula)
    fstats<-merge(data.frame(gene=geneIDs),fstats, by="gene")
    fstats<-fstats[,-grep("AveExpr", colnames(fstats))]
    return(list(stats=all.out, voom=v, fstats=fstats))
  }
}
