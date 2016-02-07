#' @title Make venn diagrams based on significance classification
#'
#'
#' @description
#' \code{counts2Venn} Take a binary significance classification matrix and
#' produce two types of venn diagrams.
#'
#' @param x A dataframe or matrix containing the binary significance calls
#' (1=significant, 0=NS)
#' @param cols A vector with the column names or numbers to use for plots.
#' Must be of length <=4.
#' @param names A vector of names for each of the venn circles.
#' @param colors A vector of colors to use for each cirlce.
#' @param type The type of venn diagram to plot. Scaled = size of circles is
#' weighted. Both = both types. Any other call gives typical venn diagrams.
#' @param legx,legy Position of legend for plot type "limma" or "both"
#' @param ... additional arguments passed to plot.
#' @details given a binary significance classification matrix, run functions to
#' produce venn diagrams.
#' if scaled, runs venneuler::venneuler venn diagrams. Otherwise, runs
#' limma::vennCounts/vennDiagram

#' @return generates a plot. Does not return anything
#'
#' @examples
#' data(kidney)
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' ### Not Run ### stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' ### Not Run ### sig<-makeBinarySig(x= stats$stats, what="Pvalue")
#' ### Not Run ### counts2Venn(x=sig, cols=c(1), names=c("treatment"),
#' ### Not Run ###    colors=c("blue"),type="limma", legx=-3.3,legy=-3)
#' @importFrom  venneuler venneuler
#' @export
counts2Venn<-function(x, cols, names, colors="black", type="both",legx=0, legy=0,...){
  if(type=="both"){
    par(mfrow=c(2,1))
  }else{
    par(mfrow=c(1,1))
  }
  if(type %in% c("scaled","both")){

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
    if(length(names)>1){
      areas<-sqrt(colSums(x[,cols])/pi)
      areas<-areas/(max(areas)*.25)
      legend("bottomleft", legend=colSums(x[,cols]), col=colors, pt.cex=areas, pch=1, bty="n",
             cex=.5, inset=c(.05,.05))
    }
  }
  par(mfrow=c(1,1))
}

