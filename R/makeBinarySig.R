#' @title Convert p-values to binary significance calls
#'
#'
#' @description
#' \code{makeBinarySig} Convert some sort of significance distribution
#' to binary signifcance calls.
#'
#' @param x A dataframe or matrix containing the significance (p/q/fdr p)
#' values to convert
#' @param alpha The threshold that defines what is and is not significance
#' @param what Character vector defining a string to look for in the column names.
#' For example, "q.value"
#' @param na.include Logical, should NAs be classified as not significant? If false,
#' all rows with any NAs are dropped.
#' @param verbose Logical, should the number of significant values be printed
#' @details Given a matrix with p/q values, this function greps for the string in
#' "what" and classifies each column as
#' either significant (1) or not (0). Esspecially useful for a dataset that has many
#' p-value or q-value columns

#' @return a dataframe with all values in the original matrix with
#' column names containing what, where the values have been
#' transformed to binary
#'
#' @examples
#' \dontrun{
#' data(kidney) # from SimSeq
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-data.frame(rep=kidney$replic, treatment=kidney$treatment)
#' stats<-pipeLIMMA(counts=counts, info=info, formula = " ~ treatment", block=NULL)
#' sig<-makeBinarySig(x= stats$stats, what="Pvalue")
#' }
#' @export
makeBinarySig<-function(x, alpha=0.05, what="q.value",
                        na.include=TRUE, verbose=TRUE){
  sig.q<-data.matrix(x[,grep(what, colnames(x))])
  if(!na.include){
    sig.q<-sig.q[complete.cases(sig.q),]
  }
  if(length(grep(what, colnames(x)))==1){
    if(na.include){
      sig.q[is.na(sig.q)]<-1
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
