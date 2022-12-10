library(phytools)

# read a tree
tree <- readTrees("stoch_map_test_tmp/tree.nexus")[[1]][[1]]

# read the samples as simmaps
samples <- readTrace("stoch_map_test_tmp/maps.log")[[1]]$simmap
simmaps <- lapply(samples, function(x) phytools::read.simmap(text = x))
class(simmaps) <- "multiPhylo"

stoch_map_df <- RevGadgets:::processStochMaps(tree, 
                                              simmap = simmaps,
                                              states = as.character(0:4))

# make colors
colors <- colFun(5)
names(colors) <- as.character(0:4)

RevGadgets:::plotStochMaps(tree = tree,
                           maps = stoch_map_df,
                           color_by = "MAP",
                           colors = colors,
                           tip_labels = F)
