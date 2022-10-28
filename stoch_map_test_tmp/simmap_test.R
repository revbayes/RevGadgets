library(phytools)

# read a tree
tree <- readTrees("stoch_map_test_tmp/tree.nexus")[[1]][[1]]

# read the samples as simmaps
samples <- readTrace("stoch_map_test_tmp/maps.log")[[1]]$simmap
simmaps <- lapply(samples, function(x) read.simmap(text = x))
class(simmaps) <- "multiPhylo"

stoch_map_df <- RevGadgets:::processStochMaps(tree, simmap = simmaps,
                                              states = as.character(0:4))

# make colors
colors <- colFun(5)
names(colors) <- states

RevGadgets:::plotStochMaps(tree = tree,
                           maps = stoch_map_df,
                           color_by = "MAP",
                           colors = colors,
                           tree_layout = "circular",
                           tip_labels = F)
