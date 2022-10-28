# run "stoch_map_test_tmp/simmaps.Rev" to generate output for this test

# read a tree
tree <- readTrees("stoch_map_test_tmp/tree.nexus")[[1]][[1]]


# process samples
stoch_map_df <- RevGadgets:::processStochMaps(tree, "stoch_map_test_tmp/maps.log", 
                                              states = as.character(0:4), 
                                              burnin = 0.1)

RevGadgets:::plotStochMaps(tree = tree,
                           maps = stoch_map_df,
                           color_by = "MAP",
                           colors = "default",
                           tree_layout = "rectangular",
                           tip_labels = F)

