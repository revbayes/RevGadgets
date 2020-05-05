#' Process ancestral states (discrete)
#'
#' [ function tags to be written ]
#'
#' @export
#'

# TO DO
# - add comments/tags to new files
# - create test script with cladogenetic events
# - break main plotting function into specialized backend functions based on `summary_statistic`
# - merge how anc_state and start_state/end_state are processed
# - develop list of plotting aesthetics to support
# - eliminate user arguments when possible (particularly those relating to plot/marker dimensions/sizes)
# - expand support for plotting arbitrary # of anc state categories

# main processing function
processAncStatesDiscrete = function(path, state_labels = NULL) {

    # read in tree
    tree <- readTrees(path)
    t <- tree[[1]][[1]]

    # process column names
    include_start_states = F
    if ("anc_state_1" %in% names(t@data)) {
        # do nothing
    } else if ("start_state_1" %in% names(t@data) && "end_state_1" %in% names(t@data)) {
        include_start_states = T
    } else {
        error("tree file does not contain expected state labels: [\'anc_state\'] or [\'start_state\' and \'end_state\']")
    }

    # add state labels
    t = .assign_state_labels(t, state_labels, include_start_states)

    # add range for pp factors
    t = .set_pp_factor_range(t, include_start_states)

    # return processed TreeIO object
    return(t)
}
