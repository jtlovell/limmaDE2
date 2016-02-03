counts2Venn<-function(x, cols, names, colors=c("darkblue","green","cyan","darkred"), type="scaled",...){

  if(type=="scaled"){
    require(venneuler, warn.conflicts = FALSE, quietly=TRUE)

    mat<-as.matrix(x[,cols])
    colnames(mat)<-names
    cs<-apply(mat,2, function(x) sum(x!=0))
    c2<-lapply(data.frame(combn(colnames(mat),m=2)),function(x) as.vector(x))
    if(length(names)>1){
      cs2<-sapply(c2, function(x) {
        temp<-mat[,x]
        sum(temp[,1]!=0 & temp[,2]!=0)
      })
      names(cs2)<-sapply(c2, function(x) paste(x, collapse="_"))
    }else{
      cs2<-""
    }
    if(length(names)>2){
      c3<-lapply(data.frame(combn(colnames(mat),m=3)),function(x) as.vector(x))
      cs3<-sapply(c3, function(x) {
        temp<-mat[,x]
        sum(temp[,1]!=0 & temp[,2]!=0 & temp[,3]!=0)
      })
      names(cs3)<-sapply(c3, function(x) paste(x, collapse="_"))
    }else{
      cs3<-""
    }
    if(length(names)>3){
      c4<-lapply(data.frame(combn(colnames(mat),m=4)),function(x) as.vector(x))
      cs4<-sapply(c4, function(x) {
        temp<-mat[,x]
        sum(temp[,1]!=0 & temp[,2]!=0 & temp[,3]!=0 & temp[,4]!=0)
      })
      names(cs4)<-sapply(c4, function(x) paste(x, collapse="_"))
    }else{
      cs4<-""
    }
    cs<-paste(paste(names(cs), cs, sep="="), collapse="  ")
    cs2<-paste(paste(names(cs2), cs2, sep="="), collapse="  ")
    cs3<-paste(paste(names(cs3), cs3, sep="="), collapse="  ")
    cs4<-paste(paste(names(cs4), cs4, sep="="), collapse="  ")
    plot(v <- venneuler(mat != 0, colors), col=colors,
         sub=paste(cs,"\n",cs2,"\n",cs3,"\n",cs4), cex.sub=.8,...)
  }else{
    vc<-vennCounts(x[,cols])
    vennDiagram(vc, names=names, circle.col=colors, cex=c(1.2,.8,.5),...)
    areas<-sqrt(colSums(x[,cols])/pi)
    areas<-areas/(max(areas)*.25)
    legend(x=0,y=0, legend=colSums(x[,cols]), col=colors, pt.cex=areas, pch=1, bty="n",
           cex=.5,xjust=0,yjust=0)
  }
}
