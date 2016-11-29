#' @title Compare semantic similarity of GO.
#'
#' @description
#' \code{distGO} A wrapper for GOSemSim::mgoSim, where similarity among
#' GO terms stored in a list are returned as a matrix
#' @param go.list A named list, where each element stores a character vector of GO
#' ids. GO terms should be sorted so that the most important are first.
#' @param n The number of GO terms to evaluate from the begining of the vector.
#' @param ont The ontology type: one of "MF", "BP", and "CC" subontologies to be
#' passed to mgoSim
#' @param organism Organism with GO database, to be passed to mgoSim
#' @param measure Measure of similarity to be passed to mgoSim
#' @param combine Measure of combining similarity to be passed to mgoSim
#' @param verbose Should updates be printed to the console?
#' @param calc.distance Should a distance matrix (default) or a similarity
#' matrix be returned.
#' @param ... Not currently in use
#'
#' @details More here soon.

#' @return a matrix of (dis)similiarities
#' @examples
#' \dontrun{
#' more here soon
#' }
#' @export
distGO<-function(go.list, n = 10,
                 ont = "BP",
                 organism = "arabidopsis",
                 measure = "Wang",
                 combine = "rcmax",
                 verbose=TRUE,
                 calc.distance = TRUE){

  l = length(go.list)
  mat<-matrix(, nrow = l, ncol = l)
  ns<-names(go.list)
  colnames(mat)<-ns
  rownames(mat)<-ns

  for(i in names(go.list)){
    if(verbose) cat("calculating similarities for",i,"\n")
    for(j in names(go.list)){
      gos1 = go.list[[i]][1:n]
      gos2 = go.list[[j]][1:n]
      out<-mgoSim(gos1, gos2, ont = ont, organism = organism, measure = measure,
                  combine = combine)
      if(calc.distance){
        out<-(1-out)
      }
      mat[i,j]<-out
    }
  }
  return(mat)
}
