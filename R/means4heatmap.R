#' @export
means4heatmap<-function(v, grps=info$Treatment){
  v<-as.matrix(v)
  mean.v<-v[,0]
  for(i in unique(grps)){
    temp<-data.frame(rowMeans(v[,grps==i]))
    colnames(temp)<-i
    mean.v<-data.frame(mean.v,temp)
  }
  return(mean.v)
}
