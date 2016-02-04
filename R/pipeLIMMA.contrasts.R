#' @title A pipeline for LIMMA via a contrast marix.
#'
#'
#' @description
#' \code{pipeLIMMA.contrasts} Run a pipeline of LIMMA functions for
#' differential gene expression using a contrast matrix. Implementation is more
#' restrictive than pipeLIMMA and forces voom w/ quality weights and geneIDs
#' to be rownames of counts matrix.
#'
#' @param counts A count matrix
#' @param design A design matrix, usually created by a call from "model.matrix"
#' @param contrast matrix A matrix of contrasts, usually created by a call
#' from limma:makeContrasts
#' @param block A string that represents an individual that was repeatedly
#' measured, if NULL, runs the analysis without a blocking / duplicate correlation factor
#' @param printSig Logical, should statistical significance for each contrast be printed?
#' @param makePlots Logical, should p-value and q-value histograms be plotted?
#' @param verbose Logical, return progress updates?
#' @param ... additional arguments passed on to lmfit, for example a vector of sample weights
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{calculate normalization factors via edgeR::calcNormFactors}
##'  \item{2. }{Run limma::voom transformation}
##'  \item{3. }{Run limma:lmFit linear modeling}
##'  \item{4. }{Run limma:contrasts.fit contrast tests}
##'  \item{5. }{Run limma::ebayes statistical modeling}
##'  \item{6. }{Ouput statistics and other data}
##' }
##'
#' @return a dataframe with the statistics for each contrast and the overall f-test
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' f<-gsub("-","_",info$treatment)
#'
#' design <- model.matrix(~0+f)
#' contrast.matrix<-makeContrasts(fNon_Tumor-fTumor, levels=design)
#' lim.contrasts<-pipeLIMMA.contrasts(counts=counts, design=design, block=info$rep,
#'        contrast.matrix=contrast.matrix)
#' @importFrom  edgeR calcNormFactors DGEList
#' @importFrom  qvalue qvalue
#' @importFrom  qdap multigsub
#' @importFrom  ggplot2 ggplot
#' @importFrom  reshape2 melt
#' @export
pipeLIMMA.contrasts<-function(counts, design, contrast.matrix,
                     block=NULL, printSig=TRUE, makePlots=TRUE, verbose=TRUE, ...){

  if(is.null(block)) {
    useBlock=FALSE
  }else{
    useBlock=TRUE
  }

  geneIDs<-rownames(counts)
  if(verbose) cat("calculating normalization factors...\n")
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)

  if(verbose) cat("running voom normalization...\n")
  v <- voomWithQualityWeights(y, design=design, plot = F)

  if(useBlock){
    if(verbose) cat("calculating duplicate correlation and fitting LIMMA w/ block...\n")
    dupcor <- duplicateCorrelation(counts,design, block=as.factor(block))
    fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block),...)
    fit2 <- contrasts.fit(fit, contrast.matrix)
  }else{
    if(verbose) cat("fitting LIMMA w/o a block effect...\n")
    fit <- lmFit(v, design=design, ...)
    fit2 <- contrasts.fit(fit, contrast.matrix)
  }
  if(verbose) cat("generating statistics...\n")
  fit2<-eBayes(fit2)
  out.aov<-data.frame(gene=geneIDs,
                      sigma=fit2$sigma,
                      s2.post=fit2$s2.post,
                      Amean=fit2$Amean,
                      F=fit2$F,
                      p.value_Ftest=fit2$F.p.value,
                      q.value_Ftest=qvalue(fit2$F.p.value)$qvalue)
  coefs<-data.frame(fit2$coef); colnames(coefs)<-paste("coef_", colnames(coefs), sep="")
  ts<-data.frame(fit2$t); colnames(ts)<-paste("t_", colnames(ts), sep="")
  lods<-data.frame(fit2$lods); colnames(lods)<-paste("lod_", colnames(lods), sep="")
  p.values<-data.frame(fit2$p.value); colnames(p.values)<-paste("p.value_", colnames(p.values), sep="")
  q.values<-apply(p.values, 2, function(x) qvalue(x)$qvalue)
  colnames(q.values)<-gsub("p.","q.",colnames(q.values))
  all.out<-data.frame(out.aov,
                      coefs,
                      ts,
                      lods,
                      p.values,
                      q.values)
  ps<-all.out[,grep("p.value", colnames(all.out))]
  qs<-all.out[,grep("q.value", colnames(all.out))]
  if(makePlots){
    if(verbose) cat("plotting results...\n")
    ps$type="p.value"
    qs$type="q.value"
    colnames(ps)<-multigsub(c("p.value_","..."),c("","."), colnames(ps))
    colnames(qs)<-multigsub(c("q.value_","..."),c("","."), colnames(qs))
    tp<-rbind(ps,qs)
    tp<-melt(tp, id.var="type")
    tp$variable<-gsub("f","",tp$variable)
    print(
      ggplot(tp, aes(x=value, col=type, fill=type))+geom_histogram(binwidth=.01)+
        facet_wrap(~variable, scales="free_y") +theme_bw()+
        theme(strip.text.x = element_text(size = 8))+
        ggtitle("p- and q-value distributions")+
        scale_x_continuous(breaks=c(0,.5,1))+
        scale_y_continuous("number of genes")
    )
  }

  colnames(qs)<-multigsub(c("q.value_","..."),c("","."), colnames(qs))
  if(printSig){
    if(verbose) cat("printing counts of significant results...\n")
    for(i in colnames(qs)[-which(colnames(qs)=="type")]) cat("n.significant genes",i,sum(qs[,i]<=0.05),"\n", sep="\t")
  }
  return(all.out)
}
