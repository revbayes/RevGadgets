#' BROKEN Reroot Phylo BROKEN
#'
#' Reroots a phylogeny given an outgroup taxon or clade
#'
#' Modifies a tree object by rerooting using a specified
#' outgroup taxon or clade. Places the root at the midpoint
#' of the branch subtending the outgroup. If the input
#' contains multiple trees, all trees will be rerooted.
#'
#'
#' @param tree (list of lists of treedata objects; no default) Name of a list of lists of
#' treedata objects, such as produced by readTrees().
#'
#' @param outgroup (character, no default) Name of the outgroup(s). Either a single taxon
#' name or a character vector of length two to specify a clade; in this case the root
#' will be placed at the midpoint of the branch subtending the two taxa's MRCA. Modified
#' from phytools::reroot().
#'
#' @return returns a list of list of treedata objects, with the trees rooted.
#'
#' @examples
#' \dontrun{
#' file <- system.file("extdata", "sub_models/primates_cytb_GTR_MAP.tre", package="RevGadgets")
#' tree <- readTrees(paths = file)
#' # root with one taxon
#' tree_rooted <- rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")
#' # root with clade, specified by two taxa
#' tree_rooted <- rerootPhylo(tree = tree,
#'                            outgroup = c("Varecia_variegata_variegata",
#'                                         "Propithecus_coquereli"))
#' }
#' @export

rerootPhylo <- function(tree, outgroup) {
  # right now this function messes with the association between the
  # data and the nodes of the tree. Must figure out how to re-associate
  # the data
  if (!is.list(tree)) stop("tree should be a list of lists of treedata objects")
  if (class(tree[[1]][[1]]) != "treedata") stop("tree should be a list of lists of treedata objects")
  if (class(outgroup) != "character") stop("outgroup should be of class character")
  if (length(outgroup) > 2) stop("outgroup should contain 1 or 2 taxa names")

  for (i in 1:length(tree)) {
    for (j in 1:length(tree[[i]])) {
      t <-  tree[[i]][[j]]@phylo
      for (k in 1:length(outgroup)) {
        if (outgroup[k] %in% t$tip.label == FALSE) {
          stop(paste0("Outgroup ", outgroup[k],
                      " not found in tree.
                      Check spelling and underscores."))
          }
        }
      if (length(outgroup) == 1) {
        num <- which(t$tip.label == outgroup)
      } else { num <- ape::getMRCA(t, outgroup) }
      midpoint <- 0.5*t$edge.length[which(t$edge[,2] == num)]
      t <- phytools::reroot(t, num, position = midpoint)
      tree[[i]][[j]]@phylo <- t
    }
  }
  return(tree)
}
