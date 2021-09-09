#' DensiTree-style plot with branch-specific data
#'
#' This function plots a distribution of trees (e.g obtained from an
#' MCMC inference) with branch-specific rates or other data.
#' The plot is similar to those produced by DensiTree, i.e all the trees are
#' overlapped with each other. The data is expected to be given per node,
#' and will be associated with the branch above its corresponding node.
#' Its values are plotted as a color gradient.
#'
#' If no consensus tree is provided, a consensus tree will be computed.
#' This should avoid too many unnecessary crossings of edges.
#' Trees should be rooted, other wise the output may not be visually pleasing.
#' The \code{jitter} parameter controls whether to shift trees so that
#' they are not exactly on top of each other.
#' If \code{amount = 0}, no jitter is applied. If \code{random = TRUE},
#' the applied jitter is calculated as \code{runif(n, -amount, amount)},
#' otherwise \code{seq(-amount, amount, length=n)}, where \code{n}
#' is the number of trees.
#'
#' @param tree_files vector of tree files in NEXUS format with data attached
#' to the branches/nodes of the tree. All trees should have the same tip
#' labels (the order can change). Either \code{tree_files} or both
#' \code{trees} and \code{data} have to be specified.
#' @param burnin fraction of samples to discard from the tree files as burn-in.
#' Default 0.1.
#' @param trees multiPhylo object or list of trees in phylo format.
#' All trees should have the same tip labels (the order can change). Either
#' \code{tree_files} or both \code{trees} and \code{data} have to be specified.
#' @param data data to be plotted on the tree - expected to be a list of
#' vectors in the same order as the trees, each vector in the order of the
#' tips and nodes of the corresponding tree
#' @param data_name Only used when reading from \code{tree_files}.
#' Name of the data to be plotted, if multiple are present.
#' @param type character string specifying the type of phylogeny.
#' Options are "cladogram" (default) or "phylogram".
#' @param consensus A tree or character vector which is used to define
#' the order of the tip labels. If NULL will be calculated from the trees.
#' @param direction a character string specifying the direction of the tree.
#' Options are "rightwards" (default), "leftwards", "upwards" and "downwards".
#' @param scaleX whether to scale trees to have identical heights.
#' Default FALSE.
#' @param width width of the tree edges.
#' @param lty line type of the tree edges.
#' @param cex a numeric value giving the factor scaling of the tip labels.
#' @param font an integer specifying the type of font for the labels: 1
#' (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param tip.color color of the tip labels.
#' @param adj a numeric specifying the justification of the text strings of
#' the tip labels: 0 (left-justification), 0.5 (centering), or 1
#' (right-justification).
#' @param srt a numeric giving how much the labels are rotated in degrees.
#' @param keep_underscores whether the underscores in tip labels should be
#' written as spaces (the default) or left as they are (if TRUE).
#' @param label_offset a numeric giving the space between the nodes and the
#' tips of the phylogeny and their corresponding labels.
#' @param scale_bar whether to add a scale bar to the plot. Default TRUE.
#' @param jitter controls whether to shift trees. a list with two arguments:
#' the amount of jitter and random or equally spaced (see details below)
#' @param color_gradient range of colors to be used for the data, in order of
#' increasing values. Defaults to red to yellow to green.
#' @param alpha transparency parameter for tree colors. If NULL will be set
#' based on the number of trees.
#' @param bias bias applied to the color gradient. See
#' \code{\link[grDevices]{colorRampPalette}} for more details.
#' @param data_intervals value intervals used for the color gradient.
#' Can be given as a vector of interval boundaries or min and max values.If
#' NULL will be set based on the data.
#' @param \dots further arguments to be passed to plot.
#'
#' @references This code is adapted from the \code{\link[phangorn]{densiTree}}
#' function by Klaus Schliep \email{klaus.schliep@@gmail.com}.
#' densiTree is inspired from the
#' \href{https://www.cs.auckland.ac.nz/~remco/DensiTree/}{DensiTree}
#' program by Remco Bouckaert.
#'
#' Remco R. Bouckaert (2010) DensiTree: making sense of sets of phylogenetic
#' trees \emph{Bioinformatics}, \bold{26 (10)}, 1372-1373.
#'
#' @return No return value, produces plot in base R
#'
#' @examples
#'
#' # generate random trees & data
#' trees <- lapply(1:5, function(x) ape::rcoal(5))
#' data <- lapply(1:5, function(x) stats::runif(9, 1, 10))
#'
#' # densiTree plot
#' densiTreeWithBranchData(trees = trees, data = data, width = 2)
#'
#' # densiTree plot with different colors
#' densiTreeWithBranchData(trees = trees, data = data,
#'                         color_gradient = c("green", "blue"), width = 2)
#'
#' @export
#' @importClassesFrom tidytree treedata


densiTreeWithBranchData <-
  function(tree_files = NULL,
           burnin = 0.1,
           trees = NULL,
           data = NULL,
           data_name = NULL,
           type = "cladogram",
           consensus = NULL,
           direction = "rightwards",
           scaleX = FALSE,
           width = 1,
           lty = 1,
           cex = .8,
           font = 3,
           tip.color = 1,
           adj = 0,
           srt = 0,
           keep_underscores = FALSE,
           label_offset = 0.01,
           scale_bar = TRUE,
           jitter = list(amount = 0, random = TRUE),
           color_gradient = c("red", "yellow", "green"),
           alpha = NULL,
           bias = 1,
           data_intervals = NULL,
           ...) {
    if ((is.null(trees) || is.null(data)) && is.null(tree_files))
      stop("Please input either trees and data or vector
           of tree files in Nexus format.")

    if (is.null(trees) || is.null(data)) {
      treedata <- readTrees(tree_files, burnin = burnin, verbose = FALSE)
      treedata <- unlist(treedata, recursive = FALSE)
      trees <- lapply(treedata, function(x)
        x@phylo)
      class(trees) <- c(class(trees), "multiPhylo")
      data <- lapply(treedata, function(x) {
        if (is.null(data_name)) {
          if (ncol(x@data) == 1)
            return(dplyr::select(x@data, names(x@data)[1]))
          else
            stop("Multiple data fields found but no data_name given")
        }
        if (!data_name %in% names(x@data))
          stop("Missing data for given data_name")
        dplyr::select(x@data, data_name)
      })
      for (i in seq_along(treedata)) {
        treedata[[i]] <-
          tidytree::treedata(phylo = treedata[[i]]@phylo, data = data[[i]])
      }
    }
    else {
      data <-
        lapply(data, function(d)
          tibble::as_tibble(as.data.frame(d)))
      treedata <- lapply(seq_along(trees), function(idx) {
        tidytree::treedata(phylo = trees[[idx]], data = data[[idx]])
      })
    }

    if (is.null(alpha))
      alpha <- max(0.01, 1 / length(trees))

    if (!is.null(data_intervals)) {
      if (length(data_intervals) == 2) {
        min.data <- min(data_intervals)
        max.data <- max(data_intervals)
        data_intervals <-
          seq(min.data, max.data, 0.1 * (max.data - min.data))
      }
    }
    else {
      # obtain max and min of data range
      min.data <- min(as.numeric(unlist(data)))
      max.data <- max(as.numeric(unlist(data)))
      data_intervals <-
        seq(min.data, max.data, 0.1 * (max.data - min.data))
    }

    ## following code adapted from phangorn::densiTree, credit Klaus Schliep

    if (is.character(consensus)) {
      consensus <- ape::stree(length(consensus), tip.label = consensus)
      consensus$edge.length <- rep(1.0, nrow(consensus$edge))
    }
    if (is.null(consensus)) {
      consensus <- ape::consensus(trees, p = .5)
    }
    if (inherits(consensus, "multiPhylo"))
      consensus <- consensus[[1]]

    type <- match.arg(type, c("phylogram", "cladogram"))
    direction <-
      match.arg(direction,
                c("rightwards", "leftwards",  "upwards",
                  "downwards"))
    horizontal <- direction %in% c("rightwards", "leftwards")

    nTip <- as.integer(length(consensus$tip.label))
    consensus <- sort_tips_phylo(consensus)
    consensus <- ape::reorder.phylo(consensus, "postorder")

    maxBT <- max(get_MRCA_heights(trees))
    if (scaleX)
      maxBT <- 1.0
    label <- rev(pretty(c(maxBT, 0)))
    maxBT <- max(label)
    xy <- ape::plotPhyloCoor(consensus, direction = direction, ...)
    yy <- xy[, 2]

    plot.new()
    tl <- which.max(nchar(consensus$tip.label))
    sw <- strwidth(consensus$tip.label[tl], cex = cex) * 1.1

    if (direction == "rightwards") {
      plot.window(xlim = c(0, 1.0 + sw), ylim = c(0, nTip + 1))
      if (scale_bar)
        axis(
          side = 1,
          at = seq(0, 1.0, length.out = length(label)),
          labels = label
        )
    }
    if (direction == "leftwards") {
      plot.window(xlim = c(0 - sw, 1.0), ylim = c(0, nTip + 1))
      if (scale_bar)
        axis(
          side = 1,
          at = seq(0, 1.0, length.out = length(label)),
          labels = rev(label)
        )
    }
    if (direction == "downwards") {
      plot.window(xlim = c(0, nTip + 1), ylim = c(0 - sw, 1.0))
      if (scale_bar)
        axis(
          side = 2,
          at = seq(0, 1.0, length.out = length(label)),
          labels = rev(label)
        )
    }
    if (direction == "upwards") {
      plot.window(xlim = c(0, nTip + 1), ylim = c(0, 1.0 + sw))
      if (scale_bar)
        axis(
          side = 2,
          at = seq(0, 1.0, length.out = length(label)),
          labels = label
        )
    }
    tip_labels <- consensus$tip.label
    if (is.expression(consensus$tip.label))
      keep_underscores <- TRUE
    if (!keep_underscores)
      tip_labels <- gsub("_", " ", tip_labels)

    add_tiplabels(
      xy,
      tip_labels,
      direction,
      adj = adj,
      font = font,
      srt = srt,
      cex = cex,
      col = tip.color,
      label_offset = label_offset
    )

    tiporder <-  1:nTip
    names(tiporder) <- consensus$tip.label

    if (jitter$amount > 0) {
      if (jitter$random)
        jit <- stats::runif(length(trees),-jitter$amount, jitter$amount)
      else
        jit <- seq(-jitter$amount, jitter$amount, length = length(trees))
    }

    for (treeindex in seq_along(treedata)) {
      tmp <- sort_tips(treedata[[treeindex]])
      phylo <- tmp@phylo

      edge_dta <- dplyr::pull(tmp@data[phylo$edge[, 2], ])
      edge_col <-
        color_gradient(edge_dta,
                       intervals = data_intervals,
                       colors = color_gradient,
                       bias = bias)

      xy <-
        ape::plotPhyloCoor(phylo,
                           tip.height = 1:nTip,
                           direction = direction,
                           ...)
      xx <- xy[, 1]
      yy <- xy[, 2]

      if (horizontal) {
        if (scaleX)
          xx <- xx / max(xx)
        else
          xx <- xx / maxBT
        if (direction == "rightwards")
          xx <- xx + (1.0 - max(xx))
        if (jitter$amount > 0)
          yy <- yy + jit[treeindex]
      }
      else {
        if (scaleX)
          yy <- yy / max(yy)
        else
          yy <- yy / maxBT
        if (direction == "upwards")
          yy <- yy + (1.0 - max(yy))
        if (jitter$amount > 0)
          xx <- xx + jit[treeindex]
      }
      e1 <- phylo$edge[, 1]
      if (type == "cladogram")
        ape::cladogram.plot(
          phylo$edge,
          xx,
          yy,
          edge.color = grDevices::adjustcolor(edge_col, alpha.f = alpha),
          edge.width = width,
          edge.lty = lty
        )
      if (type == "phylogram") {
        Ntip <- min(e1) - 1L
        Nnode <- phylo$Nnode
        ape::phylogram.plot(
          phylo$edge,
          Ntip,
          Nnode,
          xx,
          yy,
          horizontal,
          edge.color = grDevices::adjustcolor(edge_col, alpha.f = alpha),
          edge.width = width,
          edge.lty = lty
        )
      }
    }
  }
