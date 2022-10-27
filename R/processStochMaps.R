#' processStochMaps
#' 
#' @param tree (treedata object; none) Output of readTrees() function
#' containing tree.
#' @param paths (vector of character strings; no default) File path(s) to
#' stochastic map trace(s).
#' @param states (vector of character strings; no default) The character
#' states.
#' @param num_intervals (numeric; default 1001) The number of intervals
#' to divide the tree in to.
#' 
#' @export
processStochMaps <- function(tree,
                             paths,
                             states,
                             num_intervals = 1000,
                             ...) {
    
    # compute the number of states
    nstates <- length(states)
    
    # read traces
    samples <- readTrace(paths, ...)
    
    # combine multiple samples together
    if ( length(samples) > 1 ) {
        samples <- combineTraces(samples)    
    } else {
        samples <- samples[[1]]
    }
    
    # compute the number of samples
    nsamples <- nrow(samples)
    
    # get the number of branches
    num_branches <- length(tree@phylo$edge.length)
    
    # create the index map
    map <- matchNodes(tree@phylo)
    
    # get the dt
    root_age  <- max(ape::branching.times(tree@phylo))
    dt <- root_age / num_intervals
    
    # loop over branches
    dfs <- vector("list", num_branches)
    for(i in 1:num_branches) {
        
        # get the branch indexes
        R_index   <- tree@phylo$edge[i,2]
        Rev_index <- as.character(map[R_index,2])
        
        # get the time points
        this_edge_length <- tree@phylo$edge.length[i]
        these_pts        <- seq(0, this_edge_length, by = dt)
        
        # get the samples
        branch_samples <- samples[,Rev_index]
        branch_samples <- gsub("{", "", branch_samples, fixed = TRUE)
        branch_samples <- gsub("}", "", branch_samples, fixed = TRUE)
        
        # split the per event
        branch_samples <- strsplit(branch_samples, ":")
        
        # get the times per state
        branch_samples <- lapply(branch_samples, function(sample) {
            sample <- do.call(rbind, strsplit(sample, ","))
            sample_states <- sample[,1]
            sample_times  <- as.numeric(sample[,2])
            names(sample_times) <- sample_states
            return(cumsum(sample_times))
        })
        
        # get the state per interval
        branch_states_per_interval <- do.call(rbind, lapply(branch_samples, function(sample) {
            match(names(sample)[findInterval(these_pts, sample) + 1], states)
        }))
        
        # compute probability of each state per interval
        branch_prob_per_state <- apply(branch_states_per_interval, 2, tabulate, nbins = nstates) / nsamples
        rownames(branch_prob_per_state) <- states
        
        # make the df
        this_df <- data.frame(index = Rev_index, bl = this_edge_length, x0 = these_pts, x1 = c(these_pts[-1], this_edge_length))
        this_df <- cbind(this_df, t(branch_prob_per_state))
        
        # store
        dfs[[i]] <- this_df
        
    }

    # combine the branches
    dfs <- do.call(rbind, dfs)
    
    # get node instead of index 
    # node is R's standard numbering for nodes
    # index is RevBayes specific 
    nodematch <- matchNodes(tree@phylo)
    colnames(nodematch) <- c("node", "index")
    nodematch$index  <- as.character(nodematch$index)
    dfs <- dplyr::full_join(dfs, nodematch, by = "index")
    dfs$index <- NULL
    
    return(dfs)
    
}





