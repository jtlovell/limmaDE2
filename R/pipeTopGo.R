#' @title Pipeline for GO enrichment analysis
#'
#' @description
#' \code{pipeTopGo} Invokes the R package topGO to run gene ontology enrichment analysis
#' Takes a vector of gene names and a topGO custom (or not) enrichment annotation.
#'
#' @param geneID2GO A go annotation (list). Each element is the name of a gene,
#' containing a vector of GO annotations for that gene
#' @param genes.of.interest A vector of genes to look for GO enrichment
#' @param nodes4table The number of GO terms to print in the top table.
#' Defaults to the number of GO terms with fisher test p.value <= 0.05.
#' @param nodes4graph The number of GO terms to display in the plot
#' #' Defaults to the number of GO terms with fisher testp.value <=0.05.
#'

#' @return a table with the significant (or top n) GO terms and a directed network
#' of GO terms
#' @examples
#' data(geneID2GO)
#' geneNames <- names(geneID2GO)
#' myInterestingGenes <- sample(geneNames, 100)
#' test<-pipeTopGo(geneID2GO=geneID2GO, genes.of.interest=myInterestingGenes)
#' @import  topGO
#' @export
pipeTopGo<-function(geneID2GO, genes.of.interest, nodes4table=NULL, nodes4graph=NULL){
  cat("compiling annotation and generating topGO data\n")
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% genes.of.interest))
  names(geneList) <- geneNames
  GO2geneID <- inverseList(geneID2GO)
  GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)

  cat("\n----\nrunning the GO enrichment analysis by ... \n")
  cat("\t 1) the fisher classic method\n")
  capture.output(resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher"), file="/dev/null")
  cat("\t 2) the Kolmogorov-Smirnov (KS) classic method\n")
  capture.output(resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks"), file="/dev/null")

  pValue.classic <- score(resultFisher)
  n.sig<-sum(pValue.classic<=0.05)
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  hist(pValue.classic, xlab="pvalue of Fisher exact test",
       main=paste("distribution of significance\nn.sig=",n.sig), breaks=20)
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  if(is.null(nodes4table)){
    nodes4table<-n.sig
  }
  if(is.null(nodes4graph)){
    nodes4graph<-n.sig
  }
  cat("compiling results ... \n")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     classicKS = resultKS,
                     orderBy = 1, ranksOf = 1, topNodes = nodes4table)
  cat("printing graph to device ... ")
  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = nodes4graph, useInfo = 'all')
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  return(allRes)
}
