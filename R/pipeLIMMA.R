pipeLIMMA<-function(counts, info, formula, block,
                    design=NA, use.qualityWeights=TRUE,
                    tests="all",geneIDs=NA, useBlock=TRUE,
                    getTopTable=FALSE, getEbayes=TRUE,
                    contrasts=NULL, simplify=TRUE, verbose=TRUE, ...){

  require(limma, warn.conflicts = FALSE, quietly=TRUE)
  require(edgeR, warn.conflicts = FALSE, quietly=TRUE)
  require(qvalue, warn.conflicts = FALSE, quietly=TRUE)

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
  if(tests=="all"){
    tests<-attr(design, "dimnames")[[2]]
  }
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
  }
  tests.out2<-do.call(cbind, tests.out)
  all.out<-cbind(data.frame(out),tests.out2)
  colnames(all.out)<-tolower(colnames(all.out))
  return(list(stats=all.out, voom=v, lmfit=fit, countsSize=y, simpleStats=simple))
}
