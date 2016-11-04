#' @title A pipeline for DESeq2
#'
#' @description
#' \code{pipeDESeq2} Run a pipeline of DESeq2 (if installed)
#' functions for differential gene expression.
#'
#' @param counts A count matrix
#' @param info An experimental design matrix
#' @param testType The type of statistical test to run.
#' Possible options are "Wald" (Default) or "LRT".
#' The latter requires the user to specify all full
#' (formula) and reduced models to test.
#' @param formula A character string that can be
#' coerced to a formula. Specify if a contrast
#' model is not desired.
#' @param reduced If testType = "LRT", a character
#' string that can be coerced to a formula that represents
#' a sub model to formula. If multiple formulae are specified,
#' the number of formulae must match that of the formula argument.
#' All reduced formulae must be sub models of the respective
#' formula. If testType = "Wald", ignored.
#' @param geneIDs The names of genes.
#' If NA, use row names from counts matrix
#' @param verbose Logical, return progress updates?
#' @param ... additional arguments to pass to DESeq.
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{DESeq's pipeline using the specified test}
##'  \item{2. }{Extraction of results (from DESeq2::results)}
##'  \item{3. }{Renaming of columns and combining results across tests}
##' }
##'
#' @return a list with 2 elements (if simple=TRUE)
##' the statsistics generated from DESeq2::results
##'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic,
#'                  treatment=kidney$treatment)
#' stats<-pipeDESeq(counts=counts, info=info,
#'    formula = " ~ treatment")
#' stats<-pipeLIMMA(counts=counts, info=info,
#'    formula = " ~ treatment",
#'    reduced= "~ 1",
#'    testType = "LRT")
#' }
#'
#' @export
pipeDESeq2<-function(counts, info, formula=NULL,
                     reduced=NULL, testType="Wald",
                     geneIDs=NA, verbose=TRUE,...){

  if(!requireNamespace("DESeq2", quietly = TRUE)){
    stop("install the DESeq2 package to run this function\n")
  }else{
    requireNamespace("DESeq2")
  }

  se<-SummarizedExperiment(assays = data.matrix(counts),
                           colData = DataFrame(info))
  if(is.na(geneIDs)){
    geneIDs<-rownames(se)
  }
  if(length(formula)>1 & testType !="LRT"){
    stop("multiple models are only permitted for testType == LRT")
  }
  dds<- DESeqDataSet(se = se, design = as.formula(formula[1]))
  if(testType=="LRT"){
    if(length(formula)!=length(reduced)){
      stop("must have the same number of full(formula) and reduced models \n")
    }
    if(verbose){
      cat("running Likelihood ratio tests for results for:\n")
    }
    if(length(reduced)==1){
      if(verbose){
        cat("\t",formula," vs. ",reduced,"\n")
      }
      des<-DESeq(dds,test=testType, reduced= as.formula(reduced))
      resAll<-data.frame(gene=geneIDs, results(des))
    }else{
      resAll<-data.frame(gene=geneIDs)
      for(i in 1:length(reduced)){
        if(verbose){
          cat("\t",formula[i]," vs. ",reduced[i],"\n")
        }
        dds<- DESeqDataSet(se = se, design = as.formula(formula[i]))
        des<-DESeq(dds,test=testType, reduced= as.formula(reduced[i]))
        res<-data.frame(results(des))
        ns<-gsub(reduced[i],"",formula[i], fixed=T)
        ns<-gsub("*",".",ns, fixed=T)
        for(i in c("+"," ")) ns<-gsub(i,"",ns, fixed=T)

        colnames(res)<-paste(ns,colnames(res), sep="_")
        resAll<-cbind(resAll, res)
      }
    }

  }else{
    des<-DESeq(dds)
    resAll<-data.frame(gene=geneIDs)
    contrasts=sapply(resultsNames(des), list)
    if(verbose){
      cat("compiling results for contrast:\n")
    }
    for(i in 1:length(contrasts)){
      if(verbose){
        cat("\t",contrasts[i][[1]],"\n")
      }
      res <- data.frame(results(des,
                                contrast=contrasts[i],
                                pAdjustMethod = "BH"))
      colnames(res)<-paste(contrasts[i][[1]], colnames(res), sep="_")
      resAll<-cbind(resAll,res)
    }
  }
  return(resAll)
}

