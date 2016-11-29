#' @title Generate an igraph object from wgcna output
#'
#' @description
#' \code{wgcna2igraph} Function to cull and plot WGCNA networks. Requires
#' igraph and WGCNA to be installed.
#' @param net WGCNA generated network. Usually from WGCNA::blockwiseModules
#' @param datExpr The expression dataset, transposed so that genes are columns
#' and individuals are rows.
#' @param modules2plot The names (usually colors) of the WGCNA modules to plot.
#' All elements of modules2plot must be in the first element of net.
#' @param colors2plot The colors for modules2plot. Must match the length of
#' modules2plot.
#' @param kME.threshold The kME threshold to retain a node
#' @param adjacency.threshold The adjacency threshold to retain an edge.
#' @param adj.power The power used to calculate the adjacency matrix
#' @param node.size If >0, plot the nodes with a given size.
#' @param frame.color If node.size > 0, can specify the color for node outlines
#' @param node.color If node.size > 0, can specify the color for nodes
#' @param edge.alpha Numeric [0-1] specifying the transparency of the edges
#' @param verbose The position of the legend. Defaults to the top right.
#' @param returnNet Should the network be returned? If FALSE, a list of the
#' genes and original module colors that were input is returned.
#' @param ... additional arguments passed on WGCNA::adjacency
#'
#' @details More here soon.

#' @return an igraph network
#' @examples
#' \dontrun{
#' library(WGCNA)
#' library(igraph)
#' data(kidney) #' from simseq
#' counts<-kidney$counts
#' counts<-counts[sample(1:nrow(counts),1000),]
#' info<-with(kidney,
#'            data.frame(id = paste(replic, treatment, sep = "_"),
#'                       rep=replic,
#'                       Treatment=ifelse(treatment == "Tumor","tumor","cntr"),
#'                       stringsAsFactors=F))
#' colnames(counts)<-info$id
#' stats <- pipeLIMMA(counts = counts,
#'                    info = info,
#'                    block = NULL,
#'                    formula = "~ Treatment")
#'
#' datExpr.1=t(stats$voom$E)
#' pow=6
#' net.1 = blockwiseModules(datExpr.1, power = pow,
#'                          maxBlockSize = 10000, deepSplit = 2,
#'                          minModuleSize = 10,
#'                          saveTOMs = FALSE,
#'                          verbose = F)
#'
#' graph<-wgcna2igraph(net = net.1, datExpr = datExpr.1,
#'                     modules2plot = c("blue","green","turquoise","brown"),
#'                     colors2plot = c("orange","darkred","cyan","cornflowerblue"),
#'                     kME.threshold = 0.5, adjacency.threshold = 0.1,
#'                     adj.power = pow, verbose = T,
#'                     node.size = 0, frame.color = NA, node.color = NA,
#'                     edge.alpha = .5, edge.width =1)
#' plot(graph)
#' }
#'
#' @export
wgcna2igraph<-function(net, datExpr, top.n.edges = NA,
                       modules2plot = NULL, colors2plot = NULL,
                       kME.threshold = .75, adjacency.threshold = 0.1,
                       adj.power = 6, verbose = T,min.edge=2,
                       node.size = 0, frame.color = NA, node.color = NA,
                       edge.alpha = .5, edge.width =1, returnNet=TRUE,...){

  if(!returnNet & is.null(modules2plot)){
    modules2plot = unique(net[[1]])
    colors2plot = unique(net[[1]])
  }

  if(returnNet & length(colors2plot) != length(modules2plot))
    stop("colors2plot and modules2plot must have the same number of elements\n")

  if(!any(sapply(modules2plot, function(x) x %in% unique(net[[1]]))))
    stop("all modules2plot must be found in the first element of net\n")

  if(ncol(datExpr)!=length(net[[1]]))
    stop("net and datExpr must contain the same number of genes\n")

  if(!requireNamespace("WGCNA", quietly = TRUE) |
     !requireNamespace("igraph", quietly = TRUE)){
    stop("install the WGCNA and igraph packages before running\n")
  }else{
    require("igraph", quietly = TRUE)
    require("WGCNA", quietly = TRUE)
  }

  gs<-colnames(datExpr)
  cols<-net[[1]]
  names(cols)<-gs

  if(kME.threshold>0){
    if(verbose) cat("using KME values to find genes with high module membership\n")
    kme<-signedKME(datExpr = datExpr, datME = net$MEs, outputColumnName = "")
    row.names(kme)<-gs
    kmes<-sapply(1:length(gs), function(x) abs(kme[gs[x],colnames(kme) == cols[x]]))
    kmes<-data.frame(genes = gs, cols = cols, kme = kmes)
    gs<-kmes$genes[kmes$kme>=kME.threshold]
    cols<-cols[kmes$kme>=kME.threshold]
    datExpr = datExpr[,gs]
  }

  if(verbose & returnNet) cat("subsetting to modules:", paste(modules2plot, collapse = ", "),"\n")
  datExpr<-datExpr[,cols %in% modules2plot]
  cols<-cols[cols %in% modules2plot]
  gs<-gs[cols %in% modules2plot]

  if(verbose) cat("culling edges by adjacency\n")
  adj_mat<-adjacency(datExpr,power=6, ...)
  if(!is.na(top.n.edges)){
    adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
  }
  adj_mat[adj_mat > adjacency.threshold] <- 1
  adj_mat[adj_mat < adjacency.threshold] <- 0
  diag(adj_mat) <- 0
  rs<-rowSums(adj_mat)
  if(verbose) cat("removing unconnected nodes\n")
  adj_mat<-adj_mat[rs>min.edge,rs>min.edge]


  if(!returnNet){
    return(list(genes = colnames(adj_mat), cols = cols[colnames(adj_mat)]))
  }else{
    if(verbose) cat("coverting to igraph format\n")
    graph.colors = sapply(cols, function(x) colors2plot[modules2plot == x])
    net <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE,
                                       mode="upper")
    net <- simplify(net, remove.multiple = T, remove.loops = T)

    edge.start <- ends(net, es=E(net), names=T)[,1]
    edge.end <- ends(net, es=E(net), names=T)[,2]
    col.start <- graph.colors[edge.start]
    col.end <- graph.colors[edge.end]

    add.alpha <- function(col, alpha=1){
      apply(sapply(col, col2rgb)/255, 2,
            function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))
    }
    is.inMod<-col.start==col.end
    E(net)$color<-ifelse(is.inMod, add.alpha(col.start,edge.alpha),rgb(0,0,0,edge.alpha))
    E(net)$width<-edge.width
    E(net)$lty<-ifelse(E(net)$color == "#00000080" | is.na(E(net)$color), 3,1)

    V(net)$size <- node.size
    V(net)$frame.color <- frame.color
    V(net)$label <- NA
    if(node.size == 0){
      node.color = NA
    }else{
      if(is.na(node.color)){
        node.color = graph.colors
      }
    }
    V(net)$color <- node.color
    E(net)$arrow.mode <- 0
    if(verbose) cat("returning a network with",length(V(net)$size),"nodes and",length(E(net)$color),"edges\n")
    return(net)
  }
}
