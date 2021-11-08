require(colorspace)
require(ggplot2)

# modified from inset
inset.revgadgets = function (tree_view, insets, width = 0.1, height = 0.1, hjust = 0, 
    vjust = 0, x = "node", pos = 0.5) 
{
    df <- tree_view$data[as.numeric(names(insets)), ]
    
    # position subviews based on tree part
    x <- match.arg(x, c("node", "branch", "edge", "parent_shoulder"))
    if (x == "node") {
        xx <- df$x
    }
    else if (x == "parent_shoulder") {
        xx <- df$x[ match(df$parent, df$node) ]
    }
    else {
        xx <- df$branch
    }
    yy <- df$y
    xx <- xx - hjust
    yy <- yy - vjust
    
    if (length(width)==1) width = rep(width, length(insets))
    if (length(height)==1) height = rep(height, length(insets))

    # add subviews
    tree_view = tree_view + ggimage:::geom_subview(subview = insets, width = width, 
        height = height, x = xx, y = yy)

    # return treeview with subviews
    return(tree_view)
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

    # exit if no state labels provided
    if (is.null(state_labels)) {
        return(t)
    }
        
    # what is the ancestral state name tag?
    if (include_start_states) {
        state_pos_str_base = c("start_state_", "end_state_")
    } else {
        state_pos_str_base = c("anc_state_")
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
        dat[[s]]$`misc.` = rem_prob
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
    codes = c(codes, "misc.")
    return(codes)
}


################################################################################
#
# @brief Function to plot ancestral states and the associated uncertainty
#        for continuous and discrete characters.
#
#        For discrete characters 'summary_statistic="MAP"' should be used,
#        and for continuous characters 'summary_statistic="mean"'. If the
#        tree and tip labels do not fit in the screen, adjust the visible
#        area using the 'xlim_visible' argument.
#
#        If 'summary_statistic="MAP"', the maximum a posteriori ancestral
#        state will be plotted on each node. The color corresponds to the
#        character state, and the size of the circle represents the posterior
#        probability of that state. Cladogenetic models that estimate 
#        ancestral states for both the beginning and end of each branch
#        are plotted by setting "include_start_states=TRUE".
#
#        Maximum a posteriori ancestral chromosome numbers can be plotted 
#        with 'summary_statistic="MAPChromosome"'. For chromosomes the
#        color represents the posterior probability and the size of the
#        circle represents the chromosome number.
#
#        For 'summary_statistic="mean"' the color represents the size of 
#        the 95% confidence interval, and the size of the cirlce represents
#        represents the mean character state.
#
# @date Last modified: 2016-09-29
# @author Will Freyman
# @version 1.0
# @since 2016-08-31, version 1.0.0
#
# @param    tree_file               character     Path to the ancestral state tree generated by RevBayes.
# @param    summary_statistic       character   The type of summary statistic to plot.
# @param    tree_layout             character   One of 'rectangular', 'slanted', 'fan', 'circular', 'radial', or 'unrooted'.
# @param    include_start_states    logical     Plot start and end ancestral states. Used for cladogenetic models.
#
#
################################################################################
plot_ancestral_states = function(tree_file, 
                                 summary_statistic="MAP", 
                                 tree_layout="rectangular",
                                 include_start_states=FALSE, 
                                 xlim_visible=c(0, 40), 
                                 ylim_visible=NULL,
                                 tip_label_size=4, 
                                 tip_label_offset=5,
                                 tip_label_italics=FALSE,
                                 tip_node_size=2,
                                 tip_node_shape=15,
                                 node_label_size=4, 
                                 node_pp_label_size=0,
                                 node_label_nudge_x=0.1,
                                 node_pp_label_nudge_x=0.1,
                                 shoulder_label_size=3, 
                                 shoulder_label_nudge_x=-0.1, 
                                 node_pie_diameter=1.10,
                                 tip_pie_diameter=1.08,
                                 pie_nudge_x=0.0,
                                 pie_nudge_y=0.0,
                                 alpha=0.5, 
                                 node_size_range=c(6, 15), 
                                 color_low="#D55E00",
                                 color_mid="#F0E442",
                                 color_high="#009E73",
                                 show_state_legend=TRUE,
                                 show_posterior_legend=TRUE,
                                 show_tree_scale=TRUE,
                                 state_labels=NULL,
                                 state_colors=NULL,
                                 title="",
                                 fig_height=7,
                                 fig_width=7,
                                 ...) { 

    if ( (summary_statistic %in% c("MAP", "mean", "MAPChromosome", "MAPRange", "PieRange", "PieState")) == FALSE ) {
        print("Invalid summary statistic.")
        return()
    }

    # read in tree
    t = read.beast(tree_file)

    # add state labels
    #print(state_labels)
    t = assign_state_labels(t, state_labels, include_start_states)

    # add range for pp factors
    t = set_pp_factor_range(t, include_start_states)
    
    # add state colors
    use_state_colors = !is.null(state_colors)
    if (!is.null(state_colors) && !is.null(state_labels))
    {
        names(state_colors) = state_labels
    }
    
    tree = attributes(t)$phylo
    n_node = ggtree:::getNodeNum(tree)

    # remove underscores from tip labels
    attributes(t)$phylo$tip.label = gsub("_", " ", attributes(t)$phylo$tip.label)
    
    if (tip_label_italics) {
        attributes(t)$phylo$tip.label = paste("italic('", attributes(t)$phylo$tip.label, "')", sep="")
    }
    
    # add tip labels
    p = ggtree(t, layout=tree_layout, ladderize=TRUE)
    p = p + geom_tiplab(size=tip_label_size, offset=tip_label_offset, parse=tip_label_italics)
       
    if (summary_statistic == "MAPChromosome") {
        
        if (include_start_states) {
            
            if (!("start_state_1" %in% colnames(attributes(t)$data))) {
                print("Start states not found in input tree.")
                return()
            }


            # set the root's start state to NA
            attributes(t)$data$start_state_1[n_node] = NA

            # add clado daughter lineage start states on "shoulders" of tree
            # get x, y coordinates of all nodes
            x = getXcoord(tree)
            y = getYcoord(tree)
            x_anc = numeric(n_node)
            node_index = numeric(n_node)
            for (i in 1:n_node) {
                if (getParent(tree, i) != 0) {
                    # if not the root, get the x coordinate for the parent node
                    x_anc[i] = x[getParent(tree, i)]
                    node_index[i] = i
                }
            }
            shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
            p = p %<+% shoulder_data
            
            # plot the states on the "shoulders"
            p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
        
            # add ancestral states as node labels
            p = p + geom_text(aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)
        
            # show ancestral states as size / posteriors as color
            p = p + geom_nodepoint(aes(colour=as.numeric(end_state_1_pp), size=as.numeric(end_state_1)), alpha=alpha)

        } else {
        
            # add ancestral states as node labels
            p = p + geom_text(aes(label=anc_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

            # show ancestral states as size / posteriors as color
            p = p + geom_nodepoint(aes(colour=as.numeric(anc_state_1_pp), size=as.numeric(anc_state_1)), alpha=alpha)

        }

        min_low = 0.0
        max_up = 1.0
        p = p + scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=0.5)
        if (show_state_legend) {
            p = p + guides(size=guide_legend("Chromosome Number"))
        } else {
            p = p + guides(size=FALSE)
        }
        if (show_posterior_legend) {
            p = p + guides(colour=guide_legend("Posterior Probability", override.aes = list(size=8)))
        } else {
            p = p + guides(colour=FALSE)
        }

    } else if (summary_statistic == "MAPRange") {
        if (!include_start_states) {
            warning("Ignoring that include_start_states is set to FALSE")
        }
        if (!("start_state_1" %in% colnames(attributes(t)$data))) {
            print("Start states not found in input tree.")
            return()
        }

        # add ancestral states as node labels
        #p = p + geom_text(aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # set the root's start state to NA
        attributes(t)$data$start_state_1[n_node] = NA

        # add clado daughter lineage start states on "shoulders" of tree
        # get x, y coordinates of all nodes
        x = getXcoord(tree)
        y = getYcoord(tree)
        x_anc = numeric(n_node)
        node_index = numeric(n_node)
        for (i in 1:n_node) {
            if (getParent(tree, i) != 0) {
                # if not the root, get the x coordinate for the parent node
                x_anc[i] = x[getParent(tree, i)]
                node_index[i] = i
            }
        }
        shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
        p = p %<+% shoulder_data
       
        # plot the states on the "shoulders"
        p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
        p = p + geom_nodepoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)
        p = p + geom_tippoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)
        
        # show tip states as color
        #print(shoulder_data)
        #print(x_anc)
        #print(c(attributes(t)$data$start_state_1,attributes(t)$data$end_state_1))
       
        p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha) 
        
        # show ancestral states as color / posteriors as size
        p = p + geom_nodepoint(aes(colour=factor(end_state_1), size=as.numeric(end_state_1_pp)), alpha=alpha)

        if (show_state_legend) {
            p = p + guides(colour=guide_legend("Range", override.aes = list(size=8), order=1))
        } else {
            p = p + guides(colour=FALSE)
        }
        
        if (show_posterior_legend) {
            p = p + guides(size=guide_legend("Posterior probability", order=2))
        } else {
            p = p + guides(size=FALSE)
        }
        
        #return(p)

    } else if (summary_statistic == "MAP") {

        if (include_start_states) {
            print("Start states not yet implemented for MAP ancestral states.")
            return()
    
        }
        if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
            anc_data = data.frame(node=names(attributes(t)$data$end_state_1), 
                                  anc_state_1=levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
                                  anc_state_1_pp=as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
            p = p %<+% anc_data
        }

        # add ancestral states as node labels
        p = p + geom_text(aes(label=anc_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # show ancestral states as color / posteriors as size
        p = p + geom_nodepoint(aes(colour=factor(anc_state_1), size=as.numeric(anc_state_1_pp)), alpha=alpha)

        pp = as.numeric( as.vector( attributes(t)$data$anc_state_1_pp) )
        #print(pp)
        
        if (!F) {
            pp_offset_range = 2*(c(min(pp), max(pp)) - 0.5)
            nd_offset_interval = node_size_range[2] - node_size_range[1]
            nd_offset = node_size_range[1]
            node_size_range = pp_offset_range * nd_offset_interval + nd_offset
            #node_size_range[1] = node_size_range[1] * min(pp) / 0.5
            #node_size_range[2] = node_size_range[2] * max(pp)
        }

        if (node_label_size == 0) {
            p = p + geom_text(aes(label=sprintf("%.02f", as.numeric(anc_state_1_pp))), hjust="left", nudge_x=node_label_nudge_x, size=node_pp_label_size)
        }
        #p = p = scale_fill_continuous(breaks=c(0.6, 0.7, 0.8, 0.9, 1.0))
        
        # show the tip values
        p = p + geom_tippoint(aes(colour=factor(anc_state_1)), size=tip_node_size, alpha=alpha, shape=tip_node_shape)
        
        # set up the legend
        if (show_state_legend) {
            p = p + guides(colour=guide_legend("State"), order=1)        
        } else {
            p = p + guides(colour=FALSE, order=2)
        }
        if (show_posterior_legend) {
            p = p + guides(size=guide_legend("Posterior Probability"), order=3)
        } else {
            p = p + guides(size=FALSE, order=4)
        }

    } else if (summary_statistic == "mean") {
    
        if (include_start_states) {
            print("Start states not implemented for mean ancestral states.")
            return()
        }

        # add ancestral states as node labels
        p = p + geom_text(aes(label=round(mean, 2)), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

        # show the size of the 95% CI as color 
        lowers = as.numeric(levels(attributes(t)$data$lower_0.95_CI))[attributes(t)$data$lower_0.95_CI]
        uppers = as.numeric(levels(attributes(t)$data$upper_0.95_CI))[attributes(t)$data$upper_0.95_CI]
        diffs = uppers - lowers
        diffs_df = data.frame(node=names(attributes(t)$data$lower_0.95_CI), diff_vals=diffs)
        p = p %<+% diffs_df 

        min_low = min(diffs, na.rm=TRUE)
        max_up = max(diffs, na.rm=TRUE)
        mid_val = min_low + (max_up - min_low) / 2.0
        p = p + scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=mid_val)
        p = p + geom_nodepoint(aes(size=mean, colour=diff_vals), alpha=alpha)

        # show the tip values
        p = p + geom_tippoint(aes(size=mean), color="grey", alpha=alpha)

        # set up the legend
        if (show_state_legend) {
            legend_text = "Mean State"
            p = p + guides(size=guide_legend(legend_text))
        } else {
            p = p + guides(size=FALSE)
        }
        if (show_posterior_legend) {
            p = p + guides(colour=guide_legend("95% CI Width", override.aes=list(size=4)))
        } else {
            p = p + guides(colour=FALSE)
        }
    } else if (summary_statistic == "PieState") {
        if (include_start_states) {
            print("Start states not yet implemented for PieState ancestral states.")
            return()
    
        }
        
        if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
            anc_data = data.frame(node=names(attributes(t)$data$end_state_1), 
                                  anc_state_1=levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
                                  anc_state_1_pp=as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
            #p = p %<+% anc_data
        }
        
        # print tips
        p = p + geom_tippoint(aes(colour=factor(anc_state_1)), size=1e-2) 

        # plot invisible node states (for legend)
        p = p + geom_nodepoint(aes(colour=factor(anc_state_1), size=0),na.rm=TRUE, alpha=0.0)
        p = p + geom_nodepoint(aes(colour=factor(anc_state_2), size=0),na.rm=TRUE, alpha=0.0)
        p = p + geom_nodepoint(aes(colour=factor(anc_state_3), size=0),na.rm=TRUE, alpha=0.0)
        
        
        # set up the legend
        if (show_state_legend) {
            p = p + guides(colour=guide_legend("State", override.aes = list(size=5)), order=1)        
        } else {
            p = p + guides(colour=FALSE, order=2)
        }
        p = p + guides(size=FALSE)
        if (use_state_colors) {
           p = p + scale_color_manual(values=state_colors, breaks=state_labels)
        }
        
        # position legend
        p = p + theme(legend.position="left")
        
        # get anc state matrices (for pie/bar charts)
        dat_state_anc = build_state_probs(t, state_labels, include_start_states)$anc

        # make pie objects
        n_tips = length(tree$tip.label)
        n_nodes = 2 * n_tips - 1
        node_idx = (n_tips+1):n_nodes
        tip_idx = 1:n_tips
        all_idx = 1:n_nodes
        pies_anc = nodepie(dat_state_anc, cols=1:(ncol(dat_state_anc)-1), color=state_colors, alpha=alpha)
        
        # print pies
        
        # build pie diameters for tips and internal nodes
        pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )
        p_node = inset.revgadgets(tree_view=p,
                                insets=pies_anc[all_idx],
                                x="node",
                                height=pd,
                                width=pd,
                                hjust=pie_nudge_x,
                                vjust=pie_nudge_y)
        
        
        # save pdf
        # ggsave(file=paste(stree_fn,".out_state.pdf",sep=""),device="pdf",height=7,width=7)
 
        return(p_node)

    } else if (summary_statistic == "PieRange") {
     
        if (!("start_state_1" %in% colnames(attributes(t)$data))) {
            print("Start states not found in input tree.")
            return()
        }

        # set the root's start state to NA
        #attributes(t)$data$start_state_1[n_node] = NA
        
        # print tips
        p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=1e-2, alpha=alpha) 
       
        # plot invisible node states (for legend)
        p = p + geom_nodepoint(aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
        p = p + geom_nodepoint(aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
        p = p + geom_nodepoint(aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)
        
        # set up the legend
        if (show_state_legend) {
            p = p + guides(colour=guide_legend("State"), order=1)        
        } else {
            p = p + guides(colour=FALSE, order=2)
        }
        p = p + guides(size=FALSE)
        p = p + guides(colour = guide_legend(override.aes = list(size=5)))
        if (use_state_colors) {
            used_states = collect_probable_states(p)
            p = p + scale_color_manual(values=state_colors, breaks=state_labels,  name="Range", limits = used_states)
        }
        p = p + theme(legend.position="left")
        
        # # MJL: to remove later
        # break_legend = F
        # if (break_legend) {
        #     p$data$x = p$data$x + (15 - max(p$data$x))
        #     x_breaks = 0:15
        #     x_labels = rep("", 16)
        #     x_labels[ c(0,5,10,15)+1 ] = c(0,5,10,15)
        #     p = p + scale_x_continuous(breaks = x_breaks, labels = rev(x_labels), sec.axis = sec_axis(~ ., breaks = 15-c(6.15, 4.15, 2.55, 1.2), labels=c("+K","+O","+M","+H") ))
        #     p = p + theme_tree2()
        #     p = p + coord_cartesian(xlim = c(0,20), expand=TRUE)
        #     p = p + labs(x="Age (Ma)")
        #     p = add_island_times(p)
        #     p = p + theme(legend.position="left", axis.line = element_line(colour = "black"))
        #     p = p + guides(colour = guide_legend(override.aes = list(size=5), nrow=6))
        # }
            
        # get anc state matrices (for pie/bar charts)
        #print(t)
        dat_state_end = build_state_probs(t, state_labels, include_start_states)$end
        dat_state_start = build_state_probs(t, state_labels, include_start_states)$start
        
        # make pie objects
        n_tips = length(tree$tip.label)
        n_nodes = 2 * n_tips - 1
        node_idx = (n_tips+1):n_nodes
        tip_idx = 1:n_tips
        all_idx = 1:n_nodes
        
        pies_end = nodepie(dat_state_end,cols=1:(ncol(dat_state_end)-1),color=state_colors,alpha=alpha)
        pies_start = nodepie(dat_state_start,cols=1:(ncol(dat_state_start)-1),color=state_colors,alpha=alpha)
        
        pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )
        
        #n_expr = options()$expressions
        #options(expressions=n_expr * 2)
        p_node =  inset.revgadgets(tree_view=p,
                            insets=pies_end[all_idx],
                            x="node",
                            height=pd,
                            width=pd,
                            hjust=pie_nudge_x,
                            vjust=pie_nudge_y)
        
        
        p_shld = inset.revgadgets(tree_view=p_node,
                                insets=pies_start,
                                x="parent_shoulder",
                                height=node_pie_diameter*0.9,
                                width=node_pie_diameter*0.9,
                                hjust=pie_nudge_x,
                                vjust=pie_nudge_y)
        
        
        p_all = p_shld + coord_cartesian(xlim = xlim_visible, ylim=ylim_visible, expand=TRUE)
        
        return(p_shld)
    } 
    
    
    if (use_state_colors) {
        #print(state_colors)
        #print(state_labels)
        p = p + scale_color_manual(values=state_colors, breaks=as.vector(state_labels))
    }
    
    p = p + scale_radius(range = node_size_range)
    p = p + theme(legend.position="left")
    
    # show title
    p = p + ggtitle(title)
    
    # set visible area
    p = p + coord_cartesian(xlim = xlim_visible, ylim=ylim_visible, expand=TRUE)
    
    return(p)
}

