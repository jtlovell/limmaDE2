pqHists<-function(x, what.p="p.value",what.q="q.value",alpha=0.05,main="significance distribution",...){
  par(mfrow=c(1,2))
  hist(x[,grep(what.p, colnames(x))], xlab="p.value",main=paste(main, "\n","n.pvalue < ",alpha," = ",sum(x[,grep(what.p, colnames(x))]<=alpha),sep=""),...)
  hist(x[,grep(what.q, colnames(x))], xlab="q.value",main=paste(main, "\n","n.qvalue < ",alpha," = ",sum(x[,grep(what.q, colnames(x))]<=alpha),sep=""),...)
  par(mfrow=c(1,1))
}
