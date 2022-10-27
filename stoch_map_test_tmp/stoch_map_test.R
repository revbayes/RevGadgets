# run "stoch_map_test_tmp/simmaps.Rev" to generate output for this test

# read a tree
tree <- readTrees("stoch_map_test_tmp/tree.nexus")[[1]][[1]]

# process samples
stoch_map_df <- RevGadgets:::processStochMaps(tree, "stoch_map_test_tmp/maps.log", states = c("0","1"), burnin = 0.1)

#plot samples:

RevGadgets:::plotStochMaps(tree = tree,
                           maps = stoch_map_df,
                           colors = c("0" = colFun(2)[1], "1" = colFun(2)[2]))


