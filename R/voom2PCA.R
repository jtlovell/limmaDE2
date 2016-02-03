voom2PCA<-function(v, info, ids, plotit=TRUE){

  pc<-prcomp(t(v))
  prop.var<-(pc$sdev)^2 / sum(pc$sdev^2)
  dat<-data.frame(info, pc$x[,1:3])
  if(plotit){
    par(mfrow=c(2,2))
    bp<-barplot(prop.var[1:5], ylab="% Variance Explained",
                main="distribution of PCA effects", xlab="PCA Axis",
                col=c("cyan","cyan","cyan","grey","grey"))
    axis(1, at=bp[,1], labels=1:5, title)
    with(dat, plot(PC1,PC2, type="n", bty="n", main="PC1 vs. PC2"))
    with(dat, text(PC1,PC2, label=ids))
    with(dat, plot(PC1,PC3, type="n", bty="n", main="PC1 vs. PC3"))
    with(dat, text(PC1,PC3, label=ids))
    with(dat, plot(PC2,PC3, type="n", bty="n", main="PC2 vs. PC3"))
    with(dat, text(PC2,PC3, label=ids))
    par(mfrow=c(1,1))
  }
  return(dat)
}
