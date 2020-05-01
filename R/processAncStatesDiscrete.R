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



# libraries
#require(ggtree)

# main processing function
processAncStatesDiscrete = function(path,
                                    state_labels=NULL) {
    # read in tree
    t = treeio::read.beast(path)

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


# Still being developed, but this will create a probability matrix
# for all internal nodes and all sampled states. The matrix will
# be appropriate for use with the pie/bar inset function in ggtree.

#build_state_probs = function(t, state_labels, include_start_states, p_threshold = 0.01) {
#    # Generates a table that stores the states for every node in a given phylogeny
#    # States that have a posterior probability below a certain threshold at a given node,
#    # will be binned into a `other` bin (this used to be called `...`,
#    # but that name causes errors in vctrs used in tidyr)
#
#    n_states = length(state_labels)
#    n_tips = length(attributes(t)$phylo$tip.label)
#    n_node = 2 * n_tips - 1
#
#    dat = list()
#
#    if (include_start_states) {
#        state_tags = c("start","end")
#    } else {
#        state_tags = c("anc")
#    }
#
#    for (s in state_tags) {
#        dat[[s]] = data.frame( matrix(0, nrow=n_node, ncol=n_states) )
#        #dat[[s]] = cbind(node=1:n_node, dat[[s]])
#
#        for (i in 1:3)
#        {
#            m = paste(s,"_state_",i,sep="")
#            pp_str = paste(m,"_pp",sep="")
#            n_tmp = as.numeric(as.vector(attributes(t)$data$node)) # node index
#            x_tmp = as.vector(attributes(t)$data[[m]])
#            pp_tmp = as.numeric(as.vector(attributes(t)$data[[pp_str]]))
#
#            for (j in 1:length(x_tmp))
#            {
#                if (!is.na(x_tmp[j])) {
#
#                    if (pp_tmp[j] > p_threshold) {
#                        k = which(x_tmp[j]==state_labels)
#                        dat[[s]][n_tmp[j], k] = pp_tmp[j]
#                    }
#                }
#            }
#        }
#
#        # format column names
#        colnames(dat[[s]])=as.vector(unlist(state_labels))
#
#        # add probs for >3rd state under other label
#        rem_prob = c()
#        for (i in 1:nrow(dat[[s]])) {
#            rem_prob[i] = 1
#            for (j in 1:length(dat[[s]][i,])) {
#                rem_prob[i] = rem_prob[i] - dat[[s]][i,j]
#            }
#        }
#        dat[[s]]$other = rem_prob
#        dat[[s]]$node = 1:n_node
#        #print(dat[[s]][250:260,])
#    }
#
#    return(dat)
#}
#
#collect_probable_states = function(p, p_threshold=0.005)
#{
#    labels = c("end_state", "start_state")
#    index = c(1,2,3)
#
#    codes = c()
#    labels_pp = c()
#    for (l in labels) {
#        for (i in index) {
#            label_index = paste(l,"_",i,sep="")
#            label_index_pp = paste(l,"_",i,"_pp",sep="")
#            index_threshold = p$data[[ label_index_pp ]] > p_threshold
#            codes = c(codes, unique( p$data[[label_index]][ index_threshold ] ))
#        }
#    }
#    codes = unique(codes)
#    codes = c(codes, "other")
#    return(codes)
#}
