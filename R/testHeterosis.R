#' @title Bin genes by heterosis categories
#'
#' @description
#' \code{testHeterosis} Following categorization of Paschold 2012 (Genome Res.),
#' bin genes by their heterosis value. Tests are accomplished by contrasts (or similar)
#' @param contrast.stats A dataframe containing the statistical output of contrasts between
#' parents and the F1. This must contain the following columns, in this order:
##' \itemize{
##'  \item{1. }{geneID: the gene identifier}
##'  \item{2. }{p1vf1.lfc: log2fold change of parent 1 vs. the f1}
##'  \item{3. }{p2vf1.lfc: log2fold change of parent 2 vs. the f1}
##'  \item{4. }{p1vp2.lfc: log2fold change of parent 1 vs. the parent 2}
##'  \item{5. }{p1vf1.fdrp: fdr-transformed p.value of parent 1 vs. the f1}
##'  \item{6. }{p2vf1.fdrp: fdr-transformed p.value of parent 2 vs. the f1}
##'  \item{7. }{p1vp2.fdrp: fdr-transformed p.value of parent 1 vs. the parent 2}
##' }
#' @param alpha FDR threshold for determining significance
#' @param parentIDs vector of length 2: the names of P1 and P2 respectively
#' @param ... Not currently in use
#'
#' @details Categorization is as follows:
#' ##' \itemize{
##'  \item{1. }{LP<H<HP, no heterosis (additive)}
##'  \item{2. }{LP=H<HP, low parent (lp)}
##'  \item{3. }{LP<H=HP, high parent(hp)}
##'  \item{4. }{LP=HP=H, no significance (no DE)}
##'  \item{5. }{LP<HP<H, above high parent, parents are diff expressed (ahb)}
##'  \item{6. }{LP=HP<H, above high parent, parents are equal (ahb)}
##'  \item{7. }{H<LP<HP, below low parent, parents are diff expressed (blp)}
##'  \item{8. }{H<LP=HP, below low parent, parents are equal (blp)}
##'  \item{9. }{Abiguous (ambig)}
##' }
#' @return a dataframe where each row is one gene with the following columns:
##' \itemize{
##'  \item{1. }{geneID: the gene identifier}
##'  \item{2. }{highParent: p1 or p2}
##'  \item{3. }{order: the ranked order of p1, f1 and p2, separated by "."}
##'  \item{4. }{category1: the category of each gene following Paschold 2012}
##'  \item{4. }{category2: simplified category1 described in parentheses above}
##' }
##'
#' @examples
#' \dontrun{
#' test<-testHeterosis(contrast.stats=contrast.stats,
#'   alpha = 0.1,
#'   parentIDs = c("FIL2","HAL2"))
#' }
#' @export
testHeterosis<-function(contrast.stats,
                        alpha = 0.05,
                        parentIDs = c("p1","p2")){
  colnames(contrast.stats)<-c("geneID",
                              "p1vf1.lfc","p2vf1.lfc","p1vp2.lfc",
                              "p1vf1.fdrp", "p2vf1.fdrp","p1vp2.fdrp")
  geneID<-contrast.stats$geneID
  wh.hp<-with(contrast.stats, ifelse(p1vp2.lfc>0,"p1","p2"))
  names(wh.hp)<-geneID

  test.ord<-function(x){
    p1f1.up=x[1]
    p2f1.up=x[2]
    p1p2.up=x[3]
    ifelse(p1f1.up & !p2f1.up,
           "P1.F1.P2",
           ifelse(!p1f1.up & p2f1.up,
                  "P2.F1.P1",
                  ifelse(p1f1.up & p2f1.up & p1p2.up,
                         "P1.P2.F1",
                         ifelse(p1f1.up & p2f1.up & !p1p2.up,
                                "P2.P1.F1",
                                ifelse(!p1f1.up &  !p2f1.up & p1p2.up,
                                       "F1.P1.P2",
                                       ifelse(!p1f1.up &  !p2f1.up & !p1p2.up,
                                              "F1.P2.P1",NA))))))
  }

  pos<-sapply(contrast.stats[,c("p1vf1.lfc","p2vf1.lfc","p1vp2.lfc")],function(x) x>0)
  ords<-apply(pos, 1, test.ord)
  names(ords)<-geneID

  hp <- ifelse(substr(ords,1,2) == "F1",substr(ords,4,5), substr(ords,1,2))

  p1.f1.ds <- with(contrast.stats,
                   ifelse(p1vf1.fdrp>alpha,0,1)*ifelse(p1vf1.lfc>0,1,-1))
  p2.f1.ds <- with(contrast.stats,
                   ifelse(p2vf1.fdrp>alpha,0,1)*ifelse(p2vf1.lfc>0,1,-1))
  p1.p2.ds <- with(contrast.stats,
                   ifelse(p1vp2.fdrp>alpha,0,1)*ifelse(p1vp2.lfc>0,1,-1))
  hp.f1<-ifelse(hp=="P1",p1.f1.ds,p2.f1.ds)
  lp.f1<-ifelse(hp=="P2",p1.f1.ds,p2.f1.ds)
  hp.lp<-ifelse(hp=="P1",p1.p2.ds,-p1.p2.ds)
  ds<-data.frame(hp.f1,lp.f1,hp.lp, stringsAsFactors=F)
  het.cat<-ifelse(hp.f1==1 & hp.lp==1 & lp.f1==(-1),
                  "additive",
                  ifelse(lp.f1==0 & hp.f1 == 1 & hp.lp == 1,
                         "lp",
                         ifelse(lp.f1==(-1) & hp.f1 == 0 & hp.lp == 1,
                                "hp",
                                ifelse(rowSums(ds)==0,
                                       "noDE",
                                       ifelse(lp.f1==(-1) & hp.f1 == (-1) & hp.lp == 1,
                                              "ahp_parentDE",
                                              ifelse(lp.f1==(-1) & hp.f1 == (-1) & hp.lp == 0,
                                                     "ahp_noparentDE",
                                                     ifelse(lp.f1== 1 & hp.f1 == 1 & hp.lp == 1,
                                                            "blp_parentDE",
                                                            ifelse(lp.f1==1 & hp.f1 == 1 & hp.lp == 0,
                                                                   "blp_noparentDE","ambig"))))))))
  hp2<-gsub("P1",parentIDs[1],hp)
  hp2<-gsub("P2",parentIDs[2],hp2)
  ords2<-gsub("P1",parentIDs[1],ords)
  ords2<-gsub("P2",parentIDs[2],ords2)
  het.cat2<-sapply(het.cat, function(x) strsplit(x,"_")[[1]][1])
  return(data.frame(geneID, highParent = hp2, order = ords2,
                    het.category1 = het.cat, het.category2=het.cat2,
                    stringsAsFactors=F))
}
