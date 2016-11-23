#' @title Run gene ontology enrichment analyses
#'
#' @description
#' \code{pipeTopGO} Methods to simplify running limma::topGO from a table-like
#' GO annotation database.
#'
#' @param genes.of.interest A character vector representing the genes that are
#' to be tested.
#' @param GO.db The GO database in tabular format. One column must contain the
#' unique gene identifier. Gene IDs must not be replicated. Multiple GO terms
#' must be separated by comma (or similar) in a single dataframe column.
#' @param GO.db.colname The name of the column that contains the GO terms
#' @param GO.db.geneIDs The name of the GO.db column that contains the unique
#' gene identifier
#' @param GO.db.sep The character that separates GO terms.
#' @param cull2genes Specify if the background to test should be a gene set
#' other than the entire GO database
#' @param output Should the output be culled so that GO terms with P
#' values equal to 1 are not returned.
#' @details More here soon.

#' @return A tabular presentation of GO terms and the resulting statistics
#' @export
pipeTopGO<-function(genes.of.interest,
                    GO.db,
                    GO.db.colname,
                    GO.db.geneIDs,
                    GO.db.sep = ",",
                    cull2genes = NULL,
                    output = "culled"){
  if(GO.db.colname %in% colnames(GO.db))
    stop("GO.db.colname must be a column name in GO.db\n")

  if(!requireNamespace("topGO", quietly = TRUE)){
    stop("install the topGO package before running\n")
  }else{
    require("topGO", quietly = TRUE)
  }

  GO.db<-lapply(1:nrow(GO.db), function(x) strsplit(GO.db[,GO.db.colname][x],GO.db.sep)[[1]])
  names(GO.db)<-GO.db[,GO.db.geneIDs]
  nas<-sapply(GO.db, function(x) is.na(x[1]))
  GO.db<-GO.db[!nas]
  geneID2GO = GO.db
  if(!is.null(cull2genes)){
    geneID2GO<-geneID2GO[cull2genes]
  }

  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% genes.of.interest))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annotationFun = annFUN.gene2GO,
                gene2GO = geneID2GO)

  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  if(output == "culled"){
    n.non0<-length(score(resultFis)[score(resultFis)!=1])
  }else{
    n.non0<-length(score(resultFis))
  }
  allRes <- data.frame(GenTable(GOdata, resultFis, topNodes = n.non0))

  colnames(allRes)[6]<-"Pvalue"
  allRes$fdr.Pvalue<-p.adjust(allRes$Pvalue, method = "fdr")
  return(allRes)
}
