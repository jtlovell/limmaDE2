#' @title Get group means
#'
#' @description
#' \code{means4heatmap} Take normalized (or not)
#' counts data and produce a matrix that
#' can be turned into a heatmap with mean
#' values for each level.
#'
#' @param v A dataframe or matrix containing the (normalized) counts
#' @param grps The experimental factor to use to calculate groups.
#' @return a matrix with means / gene.
#'
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
