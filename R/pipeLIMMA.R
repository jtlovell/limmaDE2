#' @title A pipeline for LIMMA.
#'
#' @description
#' \code{pipeLIMMA} Run a pipeline of LIMMA functions for differential gene expression.
#'
#' @param counts A count matrix
#' @param info An experimental design matrix
#' @param formula A character string that can be coerced to a formula
#' @param block A string that represents an individual that was repeatedly measured, if NULL, runs the analysis without a blocking / duplicate correlation factor
#' @param use.qualityWeights Logical, run voom with quality weights or not?
#' @param geneIDs The names of genes. If NA, use row names from counts matrix
#' @param getTopTable Logical, return toptable statistics?
#' @param getEbayes Logical, return ebayes statistics?
#' @param simplify Logical, return a element with the F-statistics from the main model?
#' @param verbose Logical, return progress updates?
#' @params ... additional arguments, not currently in use.
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{calculate normalization factors via edgeR::calcNormFactors}
##'  \item{2. }{Run limma::voom transformation}
##'  \item{3. }{Run limma:lmFit linear modeling}
##'  \item{4. }{Run limma::ebayes statistical modeling}
##'  \item{5. }{Ouput statistics and other data}
##' }
##'
#' @return a list with 4 or 5 elements (if simple=TRUE)
##' \itemize{
##'  \item{"stats"}{: the statsistics generated from ebayes, toptable, or both}
##'  \item{"voom"}{: the voom normalized counts data}
##'  \item{"lmfit"}{: the fitted model}
##'  \item{"countsSize"}{: the normalization factors for each library}
##'  \item{"simpleStats"}{: if simplify=TRUE, a dataset with the F statistics}
##' }
##'
#' @examples
#' library(SimSeq)
#' library(limmaDE2)
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=info$rep)
#' stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
pipeLIMMA<-function(counts, info, formula, block=NULL,
                    design=NA, use.qualityWeights=TRUE,
                    geneIDs=NA, getTopTable=FALSE, getEbayes=TRUE,
                    simplify=TRUE, verbose=TRUE, ...){

  require("edgeR", quietly = TRUE, warn.conflicts = FALSE)
  require("qvalue", quietly = TRUE, warn.conflicts = FALSE)

  if(is.null(block)) {
    useBlock=FALSE
  }else{
    useBlock=TRUE
  }

  if(verbose) cat("calculating normalization factors ... \n")
  if(is.na(geneIDs)){
    geneIDs<-rownames(counts)
  }
  if(is.na(design)){
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
    if(verbose) cat("fitting linear model ... \n")
    fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block))
  }else{
    if(verbose) cat("fitting linear model ... \n")
    fit <- lmFit(v, design=design)
  }
  if(verbose) cat("processing statistics ... \n")
  fit<-eBayes(fit)
  out<-data.frame(gene=geneIDs,
                  sigma=fit$sigma,
                  s2.post=fit$s2.post,
                  Amean=fit$Amean)

  tests<-attr(design, "dimnames")[[2]]

  tests.out<-lapply(tests, function(x){
    if(getEbayes & getTopTable){
      out2<-data.frame(fit$stdev.unscaled[,x],
                       fit$coefficients[,x],
                       fit$lods[,x],
                       fit$p.value[,x],
                       qvalue(fit$p.value[,x], pi0.method="bootstrap")$qvalue)
      colnames(out2)<-paste("ebayes",x,c("stdev.unscaled","coefficients","lods","p.value","q.value"),sep="_")
      out3<-data.frame(toptable(fit, p.value=1, coef=x, number=100000))
      out3<-out3[,c("logFC","t","B")]
      colnames(out3)<-paste("tt",x,colnames(out3),sep="_")
      out2<-data.frame(out2, out3)
    }else{
      if(getTopTable){
        out2<-data.frame(toptable(fit, p.value=1, coef=x, number=100000))
        out3<-out3[,c("logFC","t","B")]
        colnames(out2)<-paste("tt", x,colnames(out2),sep="_")
      }else{
        out2<-data.frame(fit$stdev.unscaled[,x],
                        fit$coefficients[,x],
                        fit$lods[,x],
                        fit$p.value[,x],
                        if(grepl("Intercept", x)){
                          NA
                        }else{
                          qvalue(fit$p.value[,x], pi0.method="bootstrap")$qvalue
                        }
        )

        colnames(out2)<-paste("ebayes",x,c("stdev.unscaled","coefficients","lods","p.value","q.value"),sep="_")
      }
    }
    out2
  })
  if(simplify){
    fit<-fit[,-1]
    simple<-data.frame(gene=geneIDs,
                       sigma=fit$sigma,
                       s2.post=fit$s2.post,
                       Amean=fit$Amean,
                       Fstat=fit$F,
                       Fpvalue=fit$F.p.value,
                       Fqvalue=qvalue(fit$F.p.value, pi0.method="bootstrap")$qvalue)
  }else{
    simple<-NULL
  }
  tests.out2<-do.call(cbind, tests.out)
  all.out<-cbind(data.frame(out),tests.out2)
  colnames(all.out)<-tolower(colnames(all.out))
  return(list(stats=all.out, voom=v, lmfit=fit, countsSize=y, simpleStats=simple))
}
