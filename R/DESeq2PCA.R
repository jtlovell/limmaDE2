DESeq2PCA<-function(counts, info, formula, factors2Plot, factors2Extract, verbose=T, ngenes=NA){
  if(verbose) cat("running DESeq rlog transformation")
  se<-SummarizedExperiment(assays = data.matrix(counts),
                           colData = DataFrame(info))
  dds<- DESeqDataSet(se = se, design = as.formula(formula))
  rl<-rlog(dds)
  if(verbose) cat("calculating PCA, plotting and saving results")
  if(is.na(ngenes)) {
    ngenes<-nrow(counts)
  }
  plot(plotPCA(rl, intgroup=factors2Plot), ntop=ngenes)
  p<-plotPCA(rl, returnData=T, intgroup=factors2Extract)
  return(p)
}
