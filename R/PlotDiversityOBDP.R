#' Plot Diversity distribution from OBDP analysis
#' 
#' Plots the probability distribution of the number of lineages through time inferred with the 
#' Occurrence Birth Death Process
#'#'
#' @param Kt_mean (data.frame; no default) The processed data.frame previously obtained with `process_OBDP_output`.
#' @param xlab (character; "Time") The label of the x-axis.
#' @param ylab (character; "Number of lineages") The label of the y-axis.
#' @param xticks.n.breaks (numeric; 5) An integer guiding the number of major breaks. 
#' @param col.Hidden (character; "dodgerblue3") The color of the hidden lineages plot line.
#' @param col.LTT (character; "gray25") The color of the LTT plot line.
#' @param col.Total (character; "forestgreen") The color of the total lineages plot line.
#' @param col.Hidden.interval (character; "dodgerblue2") The color of the credible interval lines around the hidden lineages plot.
#' @param col.Total.interval (character; "darkolivegreen4") The color of the credible interval lines around the total lineages plot.
#' @param palette.Hidden (character; c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black")) The palette of the hidden lineages plot distribution.
#' @param palette.Total (character; c("transparent", "green4", "forestgreen", "black")) The palette of the total lineages plot distribution.
#' @param line.size (numeric; 0.7) The width of the lineage plot line.
#' @param interval.line.size (numeric; 0.5) The width of the credible interval.
#' 
#' @param show.Hidden (boolean; TRUE) Whether to show the plot for hidden lineages.
#' @param show.LTT (boolean; TRUE) Whether to show the plot for observed lineages.
#' @param show.Total (boolean; TRUE) Whether to show the plot for total lineages.
#' @param show.intervals (boolean; TRUE) Whether to show the credible intervals.
#' @param show.densities (boolean; TRUE) Whether to show the diversity densities.
#' @param show.expectations (boolean; TRUE) Whether to show the diversity expectations.
#' @param use.interpolate (boolean; TRUE) Whether to interpolate densities.
#'
#' @return A ggplot object
#' 
#' @export

plotDiversityOBDP = function( Kt_mean,
                              xlab="Time",
                              ylab="Number of lineages",
                              xticks.n.breaks = 5,
                              col.Hidden = "dodgerblue3",
                              col.LTT = "gray25",
                              col.Total = "forestgreen",
                              col.Hidden.interval = "dodgerblue2",
                              col.Total.interval = "darkolivegreen4",
                              palette.Hidden = c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                              palette.Total = c("transparent", "green4", "forestgreen", "black"),
                              line.size=0.7,
                              interval.line.size=0.5,
                              show.Hidden=TRUE,
                              show.LTT=TRUE,
                              show.Total=TRUE,
                              show.intervals=TRUE,
                              show.densities=TRUE,
                              show.expectations=TRUE,
                              use.interpolate=TRUE ){
  
  N <- length(Kt_mean)-8                                                                # Maximal number of hidden lineages
  
  ## Format Kt_mean for plotting
  Kt_mean_plot <- Kt_mean %>% tidyr::pivot_longer(-c("TimePoints", "NbObservedLin", "aggregNbHiddenLin", "aggregNbTotalLin", "NbHiddenLin0.025", "NbHiddenLin0.5", "NbHiddenLin0.975", "NbTotalLin0.025", "NbTotalLin0.5", "NbTotalLin0.975"), names_to="NbHiddenLin", values_to="ProbabilityDensity")
  Kt_mean_plot$NbHiddenLin <- as.integer(Kt_mean_plot$NbHiddenLin)

  ## Get the distribution of the total number of lineages
  Kt_mean_plot$NbTotalLin <- Kt_mean_plot$NbObservedLin + Kt_mean_plot$NbHiddenLin
  Kt_mean_plot <- Kt_mean_plot[Kt_mean_plot$ProbabilityDensity>1e-3,]                   # Remove lowest probability rows
  
  ## Plot densities and LTTs
  cols    <- c( "c1" = col.Hidden, "c2" = col.LTT, "c3" = col.Total, "H.95%" = col.Hidden.interval, "T.95%" = col.Total.interval )
  
  p <- ggplot(Kt_mean_plot, aes(x=TimePoints, y=NbTotalLin, z = ProbabilityDensity))
  
  ### Plot the number of hidden lineages : full distribution + aggregated expectation + credible interval
  if (show.Hidden){
    if (show.densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbHiddenLin, interpolate = use.interpolate,
                        fill = scales::colour_ramp(c(colorRampPalette(palette.Hidden)(N)))(Kt_mean_plot$ProbabilityDensity))
    }
    if (show.expectations){
      p <- p + geom_line(aes(y=aggregNbHiddenLin, color="c1"), size=line.size)
    }
    if (show.intervals){
      p <- p + geom_line(aes(y=NbHiddenLin, color="H.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.025, color=cols["H.95%"], linetype="twodash", size=interval.line.size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.975, color=cols["H.95%"], linetype="twodash", size=interval.line.size)
    }
  }
  
  ### Plot the total number lineages : full distribution + aggregated expectation + credible interval
  if (show.Total){
    if (show.densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbTotalLin, interpolate = use.interpolate,
                        fill = scales::colour_ramp(colorRampPalette(palette.Total)(N))(Kt_mean_plot$ProbabilityDensity))
      }
    if (show.expectations){
      p <- p + geom_line(aes(y=aggregNbTotalLin, color="c3"), size=line.size)
    }
    if (show.intervals){
      p <- p + geom_line(aes(y=NbTotalLin0.025, color="T.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.025, color=cols["T.95%"], linetype="twodash", size=interval.line.size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.975, color=cols["T.95%"], linetype="twodash", size=interval.line.size)
      }
  }
  
  ### Plot the number of observed lineages : mean
  if (show.LTT){
    p <- p +
      geom_line(aes(y=NbObservedLin, color="c2"), size=line.size)
  }
  
  ### Legend and axes
  p <- p + 
    scale_color_manual(name = "Lineages", breaks = c("c1", "H.95%", "c2", "c3", "T.95%"), 
                       values = cols, labels = c("Hidden", "95% credible interval", "LTT", "Total", "95% credible interval")) +
    scale_x_continuous(name = xlab, expand = c(0.01,0.01), n.breaks = xticks.n.breaks, minor_breaks=NULL) +   
    scale_y_continuous(name = ylab, expand = c(0.01,0.01)) + 
    theme(panel.background=element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Probability density of the number of lineages through time")
  
  return (p)
}



#' Read OBDP outputs
#' 
#' Reads and formats the outputs of an analysis with the Occurrence Birth Death Process (MCMC parameter 
#' inference + diversity estimation)
#'
#' @param start_time_trace_file               (character; no default)  Trace of the starting times along the MCMC chain.
#' @param popSize_distribution_matrices_file  (character; no default)  Kt matrices computed with `fnInferAncestralPopSize` in RevBayes.
#' @param trees_trace_file                    (character; no default)  Trace of the trees.
#' 
#' @return A data.frame
#' 
#' @export

readOBDP = function( start_time_trace_file, 
                     popSize_distribution_matrices_file,
                     trees_trace_file ){
  
  ## Import start times
  start_time_trace_lines <- readLines(start_time_trace_file)
  start_time_trace_lines <- gsub("\\[| |\\]| ,|\t|;", "", start_time_trace_lines)   # Remove unwanted characters
  start_times <- -t(read.csv(text=start_time_trace_lines[2], header=F))
    
  ## Import Kt : probability distribution of the number of hidden lineages through time
  Kt_trace_lines <- readLines(popSize_distribution_matrices_file)
  Kt_trace_lines <- gsub("\\[| |\\]| ,|\t|;", "", Kt_trace_lines)                   # Remove unwanted characters
  Kt_trace <- read.csv(text = Kt_trace_lines[-1], header = FALSE, na.strings = "nan")
  # Kt_trace[is.na(Kt_trace)] <- 0                                                  # Set NA values to 0
  S <- length(Kt_trace[,1])/length(start_times)                                     # Number of time points (nb of rows / nb of trees)
  N <- length(Kt_trace)-1                                                           # Maximal number of hidden lineages
  names(Kt_trace) <- 0:N                                                            # Set names to the number of hidden lineages
  
  ## Import the corresponding tree : get the number of observed lineages through time (LTT)
  trees_trace <- read.table(trees_trace_file, header = T)
  trees_trace$obd_tree <- lapply(trees_trace$obd_tree, function(tree){read.tree(text=as.character(tree))})
  burnin <- length(trees_trace$Iteration)-length(start_times)                       # Number of trees in the burnin
  print (paste("Burnin of", burnin, "trees over", length(trees_trace$Iteration)))
  if (burnin != 0){
    trees <- trees_trace[-(1:burnin),]                                              # Remove trees in the burnin
  }else{
    trees <- trees_trace
  }
  nb_trees <- length(trees$Iteration)                                               # Total number of trees
  iterations <- trees$Iteration
  
  Kt_trace$Iteration <- rep(iterations, each=S)                                     # Add an iteration number column
  
  ## Remove iterations with only NAs
  it_NAs <- c()
  for (it in iterations){
    if (all(is.na(Kt_trace[Kt_trace$Iteration == it,-(N+2)]))){
      print(paste("Remove iteration", it, ": contains only NAs"))
      it_NAs <- c(it_NAs, it)
    }
  }
  Kt_trace <- Kt_trace[!(Kt_trace$Iteration %in% it_NAs),]
  trees <- trees[!(trees$Iteration %in% it_NAs),]
  nb_trees <- nb_trees - length(it_NAs)
  if (!is.null(it_NAs)){start_times <- start_times[-(it_NAs+1-burnin)]}
  
  NA_rows <- which(rowSums(is.na(Kt_trace))==N+1)                                   # Get the remaining rows with only NAs
  Kt_trace[NA_rows,-(N+2)] <- cbind(rep(1, length(NA_rows)),
                                    matrix(0, nrow=length(NA_rows), ncol=N))        # These rows are considered to have 0 hidden lineages
  
  ## Add root edges lengths
  get_root_age <- function(tree){ape::ltt.plot.coords(tree)[2]}
  root_times <- sapply(trees$obd_tree, get_root_age)
  root_edge_lengths <- root_times-start_times
  for (i in 1:nb_trees){
    trees$obd_tree[[i]]$root.edge <- trees$obd_tree[[i]]$root.edge + root_edge_lengths[i]
  }
    
  ## Browse all iterations
  Kt_mean <- matrix(0, nrow=S, ncol=N+1)
  observedLin_mean <- rep(0, S)
  timePoints <- seq(0, min(start_times), length.out = S)                            # Get S time points between 0 and the oldest starting time
  for (i in 1:nb_trees){
    ### Increment the distribution of number of hidden lineages
    it <- trees$Iteration[i]
    print(paste("Tree", it, "out of", trees$Iteration[nb_trees]))
    Kt <- Kt_trace[Kt_trace$Iteration==it,-which(names(Kt_trace)=="Iteration")]
    
    if (any(Kt$`0`!=1)){                                      # Remove iterations with hidden lineages (only NA)
      Kt <- Kt/rowSums(Kt)                                    # Normalise lines to 1 (several lines at 0.9999 or 1.0001)
      Kt_mean <- Kt_mean + Kt/nb_trees
      
      ### Get the LTT coordinates
      obd_tree <- ape::collapse.singles(trees$obd_tree[[i]])       # Remove single nodes (ie. sampled ancestors)
      LTT <- data.frame(ape::ltt.plot.coords(obd_tree))            # Extract LTT coordinates
      LTT$time <- round(LTT$time, 4)                          # Reduce precision (extant tips wrongly at time -0.000001)
      
      ### Increment number of observed lineages
      getNbObservedLin <- function(t){
        if (t < LTT$time[1]) return (0)
        first_consecutive_time_point <- which(LTT$time >= t)[1]
        return (LTT$N[first_consecutive_time_point])
      }
      observedLin <- sapply(timePoints, getNbObservedLin)
      observedLin_mean <- observedLin_mean + observedLin/nb_trees
    }
    else{
      print("Error : 0 hidden lineage")
    }
  }

  ## Get the most probable number of hidden lineages (weighted mean according to their respective probabilities)
  hiddenLin <- data.frame(weightedMean=matrix(apply(Kt_mean, 1, function(x){weighted.mean(as.integer(names(Kt_mean)), x)})),
                          maxProb=matrix(apply(Kt_mean, 1, function(x){as.integer(names(Kt_mean)[which(x==max(x))])})))
  Kt_mean$aggregNbHiddenLin <- hiddenLin$weightedMean
  # Kt_mean$aggregNbHiddenLin[is.na(Kt_mean$aggregNbHiddenLin)] <- 0           # Put NA values at 0 (`weighted.mean` artifact at t0)

  ## Get the aggregated number of observed and total lineages
  Kt_mean$aggregNbTotalLin <- observedLin_mean + Kt_mean$aggregNbHiddenLin
  Kt_mean$NbObservedLin <- round(observedLin_mean)
  
  ## Get the 95% credible interval around the number of hidden lineages
  Kt_mean[1:(N+1)] <- t(apply(Kt_mean[1:(N+1)], 1, function(row){row/sum(row)}))   # Force lines to sum to 1 (correct numerical uncertainties)
  Kt_mean_cumsum <- apply(Kt_mean[1:(N+1)], 1, cumsum)
  Kt_mean$NbHiddenLin0.025 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.025)[1]-1})
  Kt_mean$NbHiddenLin0.5 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.5)[1]})
  Kt_mean$NbHiddenLin0.975 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.975)[1]})

  ## Get the 95% credible interval around the total number of lineages
  Kt_mean$NbTotalLin0.025 <- observedLin_mean + Kt_mean$NbHiddenLin0.025
  Kt_mean$NbTotalLin0.5 <- observedLin_mean + Kt_mean$NbHiddenLin0.5
  Kt_mean$NbTotalLin0.975 <- observedLin_mean + Kt_mean$NbHiddenLin0.975
  
  ## Get the aggregated number of observed and total lineages
  Kt_mean$TimePoints <- timePoints

  return(Kt_mean)
}

