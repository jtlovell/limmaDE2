#' @title Get pi0 from P-values
#'
#' @description
#' \code{reportPi0} Calculate pi0 from p-values
#' using the qvalue package
#'
#' @param ps a vector, data.frame or matrix of pvalues
#'
#' @return a vector (if is.numeric(ps)) or a
#' matrix of pi0 and the number of affected elements
#'
#' @examples
#' x <- sample(seq(0,1,by=0.0001),1000)
#' reportPi0(x)
#' df <- data.frame(x=sample(seq(0,1,by=0.0001),1000),
#'    y=sample(seq(0,.1,by=0.0001),1000),
#'    z=sample(seq(0,1,by=0.0001),1000))
#' @importFrom  qvalue qvalue
#' @export

reportPi0<-function(ps){
  if(is.numeric(ps)){
    q<-round(qvalue(ps)$pi0,2)
    n<-round((1-q)*length(ps),0)
    return(c(q,n))
  }else{
    out<-sapply(ps, function(x) {
      q<-round(qvalue(x)$pi0,2)
      n<-round((1-q)*length(x),0)
      return(c(q,n))
    })
  }
  rownames(out)<-c("pi0","n")
  return(out)
}
