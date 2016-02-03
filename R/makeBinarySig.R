makeBinarySig<-function(x, alpha=0.05, what="q.value", na.include=TRUE, verbose=TRUE){
  sig.q<-data.matrix(x[,grep(what, colnames(x))])
  if(!na.include){
    sig.q<-sig.q[complete.cases(sig.q),]
  }
  if(length(grep(what, colnames(x)))==1){
    if(na.include){
      sig.q[,i][is.na(sig.q[,i])]<-1
    }
    sig.q<-ifelse(sig.q<=alpha,1,0)
  }else{
    for(i in colnames(sig.q)) {
      if(na.include){
        sig.q[,i][is.na(sig.q[,i])]<-1
      }
      sig.q[,i]<-ifelse(sig.q[,i]<=alpha,1,0)
      if(verbose){
        cat(i, sum(sig.q[,i]),"\n", sep="\t")
      }
    }
  }
  return(sig.q)
}
