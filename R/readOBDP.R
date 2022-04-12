#' Read OBDP Outputs
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
#' @examples
#'
#' \dontrun{
#' # first run readOBDP()
#' start_time_trace_file <- 
#'      system.file("extdata", "obdp/start_time_trace.p", package="RevGadgets")
#' popSize_distribution_matrices_file <- 
#'      system.file("extdata", "obdp/Kt_trace.p", package="RevGadgets")
#' trees_trace_file <- 
#'      system.file("extdata", "obdp/mcmc_OBDP_trees.p", package="RevGadgets")
#'     
#' Kt_mean <- readOBDP( start_time_trace_file=start_time_trace_file, 
#'                      popSize_distribution_matrices_file=popSize_distribution_matrices_file, 
#'                      trees_trace_file=trees_trace_file )
#' 
#' # then get the customized ggplot object with plotDiversityOBDP()
#' p <- plotDiversityOBDP( Kt_mean,
#'                         xlab="Time (My)",
#'                         ylab="Number of lineages",
#'                         xticks_n_breaks=21,
#'                         col_Hidden="dodgerblue3",
#'                         col_LTT="gray25",
#'                         col_Total="forestgreen",
#'                         col_Hidden_interval="dodgerblue2",
#'                         col_Total_interval="darkolivegreen4",
#'                         palette_Hidden=c("transparent", "dodgerblue2", "dodgerblue3", 
#'                                          "dodgerblue4", "black"),
#'                         palette_Total=c("transparent", "green4", "forestgreen", "black"),
#'                         line_size=0.7,
#'                         interval_line_size=0.5,
#'                         show_Hidden=TRUE,
#'                         show_LTT=TRUE,
#'                         show_Total=TRUE,
#'                         show_intervals=TRUE,
#'                         show_densities=TRUE,
#'                         show_expectations=TRUE,
#'                         use_interpolate=TRUE )
#' 
#' # basic plot
#' p
#' 
#' # option: add a stratigraphic scale
#' library(deeptime)
#' library(ggplot2)
#' q <- gggeo_scale(p, dat="periods", height=unit(1.3, "line"), abbrv=F, size=4.5, neg=T)
#' r <- gggeo_scale(q, dat="epochs", height=unit(1.1, "line"), abbrv=F, size=3.5, neg=T, 
#'                     skip=c("Paleocene", "Pliocene", "Pleistocene", "Holocene"))
#' s <- gggeo_scale(r, dat="stages", height=unit(1, "line"), abbrv=T, size=2.5, neg=T)
#' s
#' }
#' 
#' @export
#' @importFrom ape read.tree collapse.singles ltt.plot.coords
#' @importFrom stats weighted.mean
#' @importFrom grDevices colorRampPalette

readOBDP = function( start_time_trace_file, 
                     popSize_distribution_matrices_file,
                     trees_trace_file ){
  
  ## Enforce argument matching
  if (!is.character(start_time_trace_file)) stop("start_time_trace_file must be a single character string")
  if (!is.character(popSize_distribution_matrices_file)) stop("popSize_distribution_matrices_file must be a single character string")
  if (!is.character(trees_trace_file)) stop("trees_trace_file must be a single character string")
  
  ## Check if trace files exist
  if (!file.exists(start_time_trace_file)) stop("start_time_trace_file does not exist")
  if (!file.exists(popSize_distribution_matrices_file)) stop("popSize_distribution_matrices_file does not exist")
  if (!file.exists(trees_trace_file)) stop("trees_trace_file does not exist")
  
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
  get_root_age <- function(tree){ltt.plot.coords(tree)[2]}
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
      obd_tree <- collapse.singles(trees$obd_tree[[i]])       # Remove single nodes (ie. sampled ancestors)
      LTT <- data.frame(ltt.plot.coords(obd_tree))            # Extract LTT coordinates
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