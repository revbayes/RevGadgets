#' Process Anc States
#'
#' Process data for ancestral states plotting
#'
#' @param path (character string; no default) File path to annotated tree.
#' @param state_labels (character vector; NULL) Vector of labels for ancestral states
#' named with the current state labels in annotated tree file (as characters).
#' @param labels_as_numbers (logical; FALSE) Should the state labels be treated as
#' integers (for example, as chromosome numbers)?
#'
#' @examples
#'
#' \dontrun{
#'
#' # standard ancestral state estimation example
#' file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
#' example <- processAncStates(file,
#'                             state_labels = c("1" = "Awesome",
#'                                              "2" = "Beautiful",
#'                                              "3" = "Cool!"))
#'
#' #chromosome evolution example
#' file <- system.file("extdata", "chromo/ChromEvol_simple_final.tree", package="RevGadgets")
#' chromo_example <- processAncStates(file, labels_as_numbers = TRUE)
#' }
#'
#' @export
#'
processAncStates <- function(path, state_labels = NULL, labels_as_numbers = FALSE) {
    #recover()
    # read in tree
    tree <- readTrees(path)
    t <- tree[[1]][[1]]

    # process column names
    include_start_states = F
    if ("anc_state_1" %in% names(t@data)) {
        # do nothing
    } else if ("start_state_1" %in% names(t@data) && "end_state_1" %in% names(t@data)) {
        include_start_states <- T
    } else {
        stop("tree file does not contain expected state labels: [\'anc_state\'] or [\'start_state\' and \'end_state\']")
    }

    # add state labels
    t <- .assign_state_labels(t, state_labels, include_start_states, labels_as_numbers)

    # add range for pp factors
    t <- .set_pp_factor_range(t, include_start_states)

    # return processed TreeIO object
    return(t)
}
