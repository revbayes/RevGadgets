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
require(ggtree)

# main processing function
processAncStatesDiscrete = function(tree_file,
                                    state_labels=NULL) {
    # read in tree
    t = read.beast(tree_file)

    # process column names
    include_start_states = F
    if ("anc_state_1" %in% names(t@data)) {
        # do nothing
    } else if ("start_state_1" %in% names(t@data) && "end_state_1" %in% names(t@data)) {
        include_start_states = T
    } else {
        error("tree_file does not contain expected state labels: [\'anc_state\'] or [\'start_state\' and \'end_state\']")
    }
    
    # add state labels
    t = assign_state_labels(t, state_labels, include_start_states)

    # add range for pp factors
    t = set_pp_factor_range(t, include_start_states)

    # return processed TreeIO object
    return(t)
}


# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getXcoord2 <- function(x, root, parent, child, len, start=0, rev=FALSE) {
    x[root] <- start
    x[-root] <- NA  ## only root is set to start, by default 0

    currentNode <- root
    direction <- 1
    if (rev == TRUE) {
        direction <- -1
    }
    while(anyNA(x)) {
        idx <- which(parent %in% currentNode)
        newNode <- child[idx]
        x[newNode] <- x[parent[idx]]+len[idx] * direction
        currentNode <- newNode
    }

    return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getXcoord <- function(tr) {
    edge <- tr$edge
    parent <- edge[,1]
    child <- edge[,2]
    root <- ggtree:::getRoot(tr)

    len <- tr$edge.length

    N <- ggtree:::getNodeNum(tr)
    x <- numeric(N)
    x <- getXcoord2(x, root, parent, child, len)
    return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getYcoord <- function(tr, step=1) {
    Ntip <- length(tr[["tip.label"]])
    N <- ggtree:::getNodeNum(tr)

    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]

    cl <- split(child, parent)
    child_list <- list()
    child_list[as.numeric(names(cl))] <- cl

    y <- numeric(N)
    tip.idx <- child[child <= Ntip]
    y[tip.idx] <- 1:Ntip * step
    y[-tip.idx] <- NA

    currentNode <- 1:Ntip
    while(anyNA(y)) {
        pNode <- unique(parent[child %in% currentNode])
        ## piping of magrittr is slower than nested function call.
        ## pipeR is fastest, may consider to use pipeR
        ##
        ## child %in% currentNode %>% which %>% parent[.] %>% unique
        ## idx <- sapply(pNode, function(i) all(child[parent == i] %in% currentNode))
        idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
        newNode <- pNode[idx]

        y[newNode] <- sapply(newNode, function(i) {
            mean(y[child_list[[i]]], na.rm=TRUE)
            ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
        })

        currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
        ## currentNode <- c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
        ## parent %in% newNode %>% child[.] %>%
        ##     `%in%`(currentNode, .) %>% `!` %>%
        ##         currentNode[.] %>% c(., newNode)
    }

    return(y)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
getParent <- function(tr, node) {
    if ( node == ggtree:::getRoot(tr) )
        return(0)
    edge <- tr[["edge"]]
    parent <- edge[,1]
    child <- edge[,2]
    res <- parent[child == node]
    if (length(res) == 0) {
        stop("cannot found parent node...")
    }
    if (length(res) > 1) {
        stop("multiple parent found...")
    }
    return(res)
}

# set custom state labels
assign_state_labels = function(t, state_labels, include_start_states, n_states=3)
{

    # what is the ancestral state name tag?
    if (include_start_states) {
        state_pos_str_base = c("start_state_", "end_state_")
    } else {
        state_pos_str_base = c("anc_state_")
    }
  
    # send error if state_labels are provided without names
    if (!is.null(state_labels) && is.null(names(state_labels))) {
        error("names(state_labels) must identify all unlabeled state names in attributes(t)$data")
    }
    
    # generate state labels if none provided
    if ( is.null(state_labels) ) {
      warning("State labels not provided by user. Will be generated automatically.")
      states <- unique(unlist(attributes(t)$data[grepl(paste0("state_","[0-9]$"),names(attributes(t)$data))]))
      states <- states[!states == "NA"]
      states <- states[order(states)]
      state_labels <- list()
      for(i in 1:length(states) ) {
        state_labels[as.character(states[i])] = LETTERS[i]
      }
      state_labels["..."] <- "..."
    }

    # create list of ancestral state name tags
    state_pos_str_to_update = c(sapply(1:n_states, function(x) { paste(state_pos_str_base,x,sep="")}))

    # overwrite state labels
    for (m in state_pos_str_to_update)
    {
        # get the states
        x_state = attributes(t)$data[[m]]
        x_state = as.vector(x_state)
        x_state_valid = which( x_state != "NA" )
        x_state_invalid = which( x_state == "NA" )
        x_state_tmp = unlist(sapply(x_state, function(z) { state_labels[ names(state_labels)==z ] }))
        x_state[x_state_valid] = x_state_tmp
        x_state[x_state_invalid] = NA
        attributes(t)$data[[m]] = x_state
    }
    
    # Just add the state_labels here
    attributes(t)$state_labels <- state_labels
    
    return(t)
}

# set prob factors
set_pp_factor_range = function(t, include_start_states, n_states=1)
{

    # what is the ancestral state name tag?
    if (include_start_states) {
        state_pos_str_base = c("start_state_", "end_state_")
    } else {
        state_pos_str_base = c("anc_state_")
    }

    # create list of ancestral state name tags
    state_pos_str_to_update = c(sapply(1:n_states, function(x) { paste(state_pos_str_base,x,"_pp",sep="")}))

    # overwrite state labels
    for (m in state_pos_str_to_update)
    {
        x_state = attributes(t)$data[[m]]
        #levels(x_state) = c(levels(x_state))
        attributes(t)$data[[m]] = x_state
    }
    return(t)
}

# Still being developed, but this will create a probability matrix
# for all internal nodes and all sampled states. The matrix will
# be appropriate for use with the pie/bar inset function in ggtree.

build_state_probs = function(t, state_labels, include_start_states, p_threshold = 0.01) {
    # Generates a table that stores the states for every node in a given phylogeny
    # States that have a posterior probability below a certain threshold at a given node, will be binned into a `...` bin
    
    n_states = length(state_labels)
    n_tips = length(attributes(t)$phylo$tip.label)
    n_node = 2 * n_tips - 1

    dat = list()

    if (include_start_states) {
        state_tags = c("start","end")
    } else {
        state_tags = c("anc")
    }

    for (s in state_tags) {
        dat[[s]] = data.frame( matrix(0, nrow=n_node, ncol=n_states) )
        #dat[[s]] = cbind(node=1:n_node, dat[[s]])

        for (i in 1:3)
        {
            m = paste(s,"_state_",i,sep="")
            pp_str = paste(m,"_pp",sep="")
            n_tmp = as.numeric(as.vector(attributes(t)$data$node)) # node index
            x_tmp = as.vector(attributes(t)$data[[m]])
            pp_tmp = as.numeric(as.vector(attributes(t)$data[[pp_str]]))

            for (j in 1:length(x_tmp))
            {
                if (!is.na(x_tmp[j])) {

                    if (pp_tmp[j] > p_threshold) {
                        k = which(x_tmp[j]==state_labels)
                        dat[[s]][n_tmp[j], k] = pp_tmp[j]
                    }
                }
            }
        }

        # format column names
        colnames(dat[[s]])=as.vector(unlist(state_labels))

        # add probs for >3rd state under ... label
        rem_prob = c()
        for (i in 1:nrow(dat[[s]])) {
            rem_prob[i] = 1
            for (j in 1:length(dat[[s]][i,])) {
                rem_prob[i] = rem_prob[i] - dat[[s]][i,j]
            }
        }
        dat[[s]]$`...` = rem_prob
        dat[[s]]$node = 1:n_node
        #print(dat[[s]][250:260,])
    }

    return(dat)
}

collect_probable_states = function(p, p_threshold=0.005)
{
    labels = c("end_state", "start_state")
    index = c(1,2,3)

    codes = c()
    labels_pp = c()
    for (l in labels) {
        for (i in index) {
            label_index = paste(l,"_",i,sep="")
            label_index_pp = paste(l,"_",i,"_pp",sep="")
            index_threshold = p$data[[ label_index_pp ]] > p_threshold
            codes = c(codes, unique( p$data[[label_index]][ index_threshold ] ))
        }
    }
    codes = unique(codes)
    codes = c(codes, "...")
    return(codes)
}


make_states = function(label_fn, color_fn) {

    # generate colors for ranges
    range_color_list = read.csv(color_fn, header=T, sep=",", colClasses="character")
    
    # get area names
    area_names = unlist(sapply(range_color_list$range, function(y) { if (nchar(y)==1) { return(y) } }))

    # get state labels
    state_descriptions = read.csv(label_fn, header=T, sep=",", colClasses="character")
    
    # map presence-absence ranges to area names
    range_labels = sapply(state_descriptions$range[2:nrow(state_descriptions)],
        function(x) {
            present = as.vector(gregexpr(pattern="1", x)[[1]])
            paste( area_names[present], collapse="")
        })

    # map labels to colors 
    range_colors = range_color_list$color[ match(range_labels, range_color_list$range) ]
    
    # generate state/color labels
    idx = 1
    st_lbl = list()
    st_colors = c()
    for (j in 1:(nrow(state_descriptions)-1)) {
        st_lbl[[ as.character(j) ]] = range_labels[j]
        st_colors[j] = range_colors[j]
    }
    st_colors[ length(st_colors)+1 ] = "lightgray"
    st_lbl[["..."]] = "..."   
    
    return( list(state_labels=st_lbl, state_colors=st_colors) )
}