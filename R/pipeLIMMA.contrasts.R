anovaLIMMA<-function(counts, design, contrast.matrix,
                     block=NULL, useBlock=F, printSig=TRUE, makePlots=TRUE, verbose=T){

  require(limma, warn.conflicts = FALSE, quietly=TRUE)
  require(edgeR, warn.conflicts = FALSE, quietly=TRUE)
  require(qvalue, warn.conflicts = FALSE, quietly=TRUE)
  require(ggplot2, warn.conflicts = FALSE, quietly=TRUE)
  require(qdap, warn.conflicts = FALSE, quietly=TRUE)

  geneIDs<-rownames(counts)
  if(verbose) cat("calculating normalization factors...\n")
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)

  if(verbose) cat("running voom normalization...\n")
  v <- voomWithQualityWeights(y, design=design, plot = F)

  if(useBlock){
    if(verbose) cat("calculating duplicate correlation and fitting LIMMA w/ block...\n")
    dupcor <- duplicateCorrelation(counts,design, block=as.factor(block))
    fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block))
    fit2 <- contrasts.fit(fit, contrast.matrix)
  }else{
    if(verbose) cat("fitting LIMMA w/o a block effect...\n")
    fit <- lmFit(v, design=design)
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
    for(i in colnames(qs)) cat("n.significant genes",i,sum(qs[,i]<=0.05),"\n", sep="\t")
  }
  return(all.out)
}
