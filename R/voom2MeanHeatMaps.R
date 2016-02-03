voom2MeanHeatMaps<-function(v,grps=info$Treatment,rowids=info$ID,thresh=7, calcMeans=T, newIDs=NA, allHMtoo=T){
  opar<-par()
  require(gdata, warn.conflicts = FALSE, quietly=TRUE)
  require(dendroextras, warn.conflicts = FALSE, quietly=TRUE)
  require(Heatplus, warn.conflicts = FALSE, quietly=TRUE)

  if(calcMeans){
    v<-as.matrix(v)
    colnames(v)[startsWith(colnames(v),"X")]<-gsub("X","",colnames(v)[startsWith(colnames(v),"X")])
    mean.v<-v[,0]
    for(i in unique(grps)){
      temp<-data.frame(rowMeans(v[,as.character(rowids[grps==i])]))
      colnames(temp)<-i
      mean.v<-data.frame(mean.v,temp)
    }
  }else{
    mean.v<-v
    colnames(mean.v)<-newIDs
    mean.v<-mean.v[,order(newIDs)]
  }
  scl<-t(mean.v)
  scl2<-scale(scl, center=T, scale=F)
  std.mean.v<-t(scl2)
  hc<-hclust(dist(std.mean.v))
  d.out=colour_clusters(hc,h=thresh)
  cut.in<-cutree(hc,h=thresh)
  par(mfrow=c(1,1))
  pal<-colorRampPalette(c('dark blue','white','dark red'))
  plot(d.out)
  abline(h=thresh, lty=2)
  if(allHMtoo){
    v2<-data.matrix(std.mean.v)
    rownames(v2)<-NULL
    hm1<-annHeatmap2(x=v2, labels=list(row=F), scale="none",col=pal, legend=T)
    plot(hm1, cex=.5)
  }
  for(i in names(table(cut.in))[table(cut.in)>1]){
    temp<-std.mean.v[cut.in==i,]
    rownames(temp)<-NULL
    hm1<-annHeatmap2(x=data.matrix(temp), labels=list(row=F), scale="none",col=pal, legend=T)
    plot(hm1, cex.axis=.5, main = paste("cluster",i), col.main=rainbow(length(names(table(cut.in))[table(cut.in)>1]))[i])
  }
  return(mean.v)
  par(opar)
}
