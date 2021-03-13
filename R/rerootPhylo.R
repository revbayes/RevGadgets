#' Reroot Phylo
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

  if (!is.list(tree))
    stop("tree should be a list of lists of treedata objects")
  if (class(tree[[1]][[1]]) != "treedata")
    stop("tree should be a list of lists of treedata objects")
  if (class(outgroup) != "character")
    stop("outgroup should be of class character")
  if (length(outgroup) > 2)
    stop("outgroup should contain 1 or 2 taxa names")

  for (i in 1:length(tree)) {
    for (j in 1:length(tree[[i]])) {

      # Make tips names for tree and add to data object
      node_name <- .makeNodeNames(tree = tree[[i]][[j]]@phylo)
      tree[[i]][[j]]@data <- tree[[i]][[j]]@data[order(as.numeric(tree[[i]][[j]]@data$node)), ]
      tree[[i]][[j]]@data$node_name <- node_name$node_names

      # Check that outgroups are in tree
      for (k in 1:length(outgroup)) {
        if (outgroup[k] %in% tree[[i]][[j]]@phylo$tip.label == FALSE) {
          stop(paste0("Outgroup ", outgroup[k],
              " not found in tree. Check spelling and underscores."))
        }
      }

      # Reroot the tree
      tmp <- tree[[i]][[j]]@phylo
      t_rooted <- tmp
      if (length(outgroup) == 1) {
        num <- which(t_rooted$tip.label == outgroup)
      } else {
        num <- ape::getMRCA(t_rooted, outgroup)
      }

      midpoint <- 0.5 * t_rooted$edge.length[which(t_rooted$edge[, 2] == num)]
      t_rooted <- phytools::reroot(t_rooted, num, position = midpoint)

      # # re-index the nodes
      # old_labels <- unique(t_rooted$edge[t_rooted$edge[,2] > length(t_rooted$tip.label),2])
      # new_labels <- sort(old_labels)
      # old_edge <- as.vector(t_rooted$edge)
      # old_edge[old_edge %in% old_labels] <- new_labels[na.omit(match(old_edge, old_labels))]
      # t_rooted$edge <- matrix(old_edge, ncol=2)
      #
      # par(mfrow=c(1,2), mar=c(0,0,0,0))
      # plot(tmp, no.margin=TRUE, cex=0.5)
      # nodelabels()
      # plot(t_rooted, no.margin=TRUE, cex=0.5)
      # nodelabels()

      # Get node names for new, rerooted tree
      node_names_new <-
        data.frame(
          index     = as.character(length(t_rooted$tip.label) + t_rooted$Nnode),
          node      = as.character(length(t_rooted$tip.label) + t_rooted$Nnode),
          node_name = utils::tail(.makeNodeNames(tree = t_rooted)$node_names, 1)
        )

      tree[[i]][[j]]@data <- dplyr::full_join(tree[[i]][[j]]@data, node_names_new, by = c("index","node", "node_name"))

      # node_names_new <-
      #   data.frame(
      #     node_name = .makeNodeNames(tree = t_rooted)$node_names,
      #     node_name_op = .makeNodeNames(tree = t_rooted)$node_names_op,
      #     node_new = 1:(length(t_rooted$tip.label) + t_rooted$Nnode)
      #   )
      #
      # # Combine new node names with data to associate new tree
      # tree[[i]][[j]]@data <- dplyr::full_join(tree[[i]][[j]]@data, node_names_new, by = "node_name")
      #
      # # Some nodes now have no info
      # # Assign them the info from nodes that contain every BUT those taxa
      # for (k in 1:nrow(tree[[i]][[j]]@data)) {
      #   if (is.na(tree[[i]][[j]]@data$node_new[k])) {
      #     n <-
      #       which(tree[[i]][[j]]@data$node_name_op == tree[[i]][[j]]@data$node_name[k])
      #     tree[[i]][[j]]@data$node_new[k] <-
      #       tree[[i]][[j]]@data$node_new[n]
      #   }
      # }
      #
      # # replace old node ID with new
      # tree[[i]][[j]]@data$node <- NULL
      # tree[[i]][[j]]@data$node <- tree[[i]][[j]]@data$node_new

      # replace old tree with new
      tree[[i]][[j]]@phylo <- t_rooted

    }

  }

  return(tree)
}
