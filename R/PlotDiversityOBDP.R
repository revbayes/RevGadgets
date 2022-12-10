#' Plot Diversity Distribution from OBDP Analysis
#' 
#' Plots the probability distribution of the number of lineages through time inferred with the 
#' Occurrence Birth Death Process
#'#'
#' @param Kt_mean (data.frame; no default) The processed data.frame (output of readOBDP()).
#' @param xlab (character; "Time") The label of the x-axis.
#' @param ylab (character; "Number of lineages") The label of the y-axis.
#' @param xticks_n_breaks (numeric; 5) An integer guiding the number of major breaks. 
#' @param col_Hidden (character; "dodgerblue3") The color of the hidden lineages plot line.
#' @param col_LTT (character; "gray25") The color of the LTT plot line.
#' @param col_Total (character; "forestgreen") The color of the total lineages plot line.
#' @param col_Hidden_interval (character; "dodgerblue2") The color of the credible interval lines around the hidden lineages plot.
#' @param col_Total_interval (character; "darkolivegreen4") The color of the credible interval lines around the total lineages plot.
#' @param palette_Hidden (character; c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black")) The palette of the hidden lineages plot distribution.
#' @param palette_Total (character; c("transparent", "green4", "forestgreen", "black")) The palette of the total lineages plot distribution.
#' @param line_size (numeric; 0.7) The width of the lineage plot line.
#' @param interval_line_size (numeric; 0.5) The width of the credible interval.
#' 
#' @param show_Hidden (boolean; TRUE) Whether to show the plot for hidden lineages.
#' @param show_LTT (boolean; TRUE) Whether to show the plot for observed lineages.
#' @param show_Total (boolean; TRUE) Whether to show the plot for total lineages.
#' @param show_intervals (boolean; TRUE) Whether to show the credible intervals.
#' @param show_densities (boolean; TRUE) Whether to show the diversity densities.
#' @param show_expectations (boolean; TRUE) Whether to show the diversity expectations.
#' @param use_interpolate (boolean; TRUE) Whether to interpolate densities.
#'
#' @return A ggplot object
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
#' @importFrom scales colour_ramp
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes aes_ geom_line annotate scale_color_manual scale_x_continuous scale_y_continuous theme element_line element_rect element_blank ggtitle

plotDiversityOBDP = function( Kt_mean,
                              xlab="Time",
                              ylab="Number of lineages",
                              xticks_n_breaks = 5,
                              col_Hidden = "dodgerblue3",
                              col_LTT = "gray25",
                              col_Total = "forestgreen",
                              col_Hidden_interval = "dodgerblue2",
                              col_Total_interval = "darkolivegreen4",
                              palette_Hidden = c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                              palette_Total = c("transparent", "green4", "forestgreen", "black"),
                              line_size=0.7,
                              interval_line_size=0.5,
                              show_Hidden=TRUE,
                              show_LTT=TRUE,
                              show_Total=TRUE,
                              show_intervals=TRUE,
                              show_densities=TRUE,
                              show_expectations=TRUE,
                              use_interpolate=TRUE ){
  
  ## Enforce argument matching
  if (!is.data.frame(Kt_mean)) stop("Kt_mean must be a data.frame (obtained with readOBDP())")
  if (!is.character(xlab)) stop("xlab must be a single character string")
  if (!is.character(ylab)) stop("ylab must be a single character string")
  if (!is.numeric(xticks_n_breaks)) stop("xticks_n_breaks must be a numeric value")
  if (!is.character(col_Hidden)) stop("col_Hidden must be a single character string")
  if (!is.character(col_LTT)) stop("col_LTT must be a single character string")
  if (!is.character(col_Total)) stop("col_Total must be a single character string")
  if (!is.character(col_Hidden_interval)) stop("col_Hidden_interval must be a single character string")
  if (!is.character(col_Total_interval)) stop("col_Total_interval must be a single character string")
  if (!is.character(palette_Hidden)) stop("palette_Hidden must be a vector of strings")
  if (!is.character(palette_Total)) stop("palette_Total must be a vector of strings")
  if (!is.numeric(line_size)) stop("line_size must be a numeric value")
  if (!is.numeric(interval_line_size)) stop("interval_line_size must be a numeric value")
  if (!is.logical(show_Hidden)) stop("show_Hidden must be a boolean")
  if (!is.logical(show_LTT)) stop("show_LTT must be a boolean")
  if (!is.logical(show_Total)) stop("show_Total must be a boolean")
  if (!is.logical(show_intervals)) stop("show_intervals must be a boolean")
  if (!is.logical(show_densities)) stop("show_densities must be a boolean")
  if (!is.logical(show_expectations)) stop("show_expectations must be a boolean")
  if (!is.logical(use_interpolate)) stop("use_interpolate must be a boolean")
  
  ## Check if KT_mean is formatted correctly
  if (!all(c("aggregNbHiddenLin", "aggregNbTotalLin", "NbObservedLin", "NbHiddenLin0.025", 
             "NbHiddenLin0.5", "NbHiddenLin0.975", "NbTotalLin0.025", "NbTotalLin0.5", 
             "NbTotalLin0.975", "TimePoints") %in% names(Kt_mean))) stop("Kt_mean have been formatted obtained with readOBDP())")
  
  ## Format Kt_mean for plotting
  Kt_mean_plot <- pivot_longer(Kt_mean, -c("TimePoints", "NbObservedLin", "aggregNbHiddenLin", "aggregNbTotalLin", "NbHiddenLin0.025", "NbHiddenLin0.5", "NbHiddenLin0.975", "NbTotalLin0.025", "NbTotalLin0.5", "NbTotalLin0.975"), names_to="NbHiddenLin", values_to="ProbabilityDensity")
  Kt_mean_plot$NbHiddenLin <- as.integer(Kt_mean_plot$NbHiddenLin)

  ## Get the distribution of the total number of lineages
  Kt_mean_plot$NbTotalLin <- Kt_mean_plot$NbObservedLin + Kt_mean_plot$NbHiddenLin
  Kt_mean_plot <- Kt_mean_plot[Kt_mean_plot$ProbabilityDensity>1e-3,]                   # remove lowest probability rows
  
  ## Plot densities and LTTs
  cols    <- c( "c1" = col_Hidden, "c2" = col_LTT, "c3" = col_Total, "H.95%" = col_Hidden_interval, "T.95%" = col_Total_interval )
  
  p <- ggplot(Kt_mean_plot, aes_(x=~TimePoints, y=~NbTotalLin, z=~ProbabilityDensity))
  
  ### Plot the number of hidden lineages : full distribution + aggregated expectation + credible interval
  N <- length(Kt_mean)-8       # maximum number of hidden lineages
  if (show_Hidden){
    if (show_densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbHiddenLin, interpolate = use_interpolate,
                        fill = colour_ramp(c(colorRampPalette(palette_Hidden)(N)))(Kt_mean_plot$ProbabilityDensity))
    }
    if (show_expectations){
      p <- p + geom_line(aes_(y=~aggregNbHiddenLin, color="c1"), linewidth=line_size)
    }
    if (show_intervals){
      p <- p + geom_line(aes_(y=~NbHiddenLin, color="H.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.025, color=cols["H.95%"], linetype="twodash", linewidth=interval_line_size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.975, color=cols["H.95%"], linetype="twodash", linewidth=interval_line_size)
    }
  }
  
  ### Plot the total number lineages : full distribution + aggregated expectation + credible interval
  if (show_Total){
    if (show_densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbTotalLin, interpolate = use_interpolate,
                        fill = colour_ramp(colorRampPalette(palette_Total)(N))(Kt_mean_plot$ProbabilityDensity))
      }
    if (show_expectations){
      p <- p + geom_line(aes_(y=~aggregNbTotalLin, color="c3"), linewidth=line_size)
    }
    if (show_intervals){
      p <- p + geom_line(aes_(y=~NbTotalLin0.025, color="T.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.025, color=cols["T.95%"], linetype="twodash", linewidth=interval_line_size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.975, color=cols["T.95%"], linetype="twodash", linewidth=interval_line_size)
      }
  }
  
  ### Plot the number of observed lineages : mean
  if (show_LTT){
    p <- p + geom_line(aes_(y=~NbObservedLin, color="c2"), linewidth=line_size)
  }
  
  ### Legend
  breaks <- c()
  labs <- c()
  if (show_Hidden){
    breaks <- c(breaks, "c1"); labs <- c(labs, "Hidden")
    if (show_intervals) {breaks <- c(breaks, "H.95%"); labs <- c(labs, "95% credible interval")}
  }
  if (show_LTT) {breaks <- c(breaks, "c2"); labs <- c(labs, "LTT")}
  if (show_Total){
    breaks <- c(breaks, "c3"); labs <- c(labs, "Total")
    if (show_intervals) {breaks <- c(breaks, "T.95%"); labs <- c(labs, "95% credible interval")}
  }
  p <- p + scale_color_manual(name = "Lineages", breaks = breaks, values = cols, labels = labs)
  
  ### Axes
  p <- p + 
    scale_x_continuous(name = xlab, expand = c(0.01,0.01), n.breaks = xticks_n_breaks, minor_breaks=NULL) +   
    scale_y_continuous(name = ylab, expand = c(0.01,0.01)) + 
    theme(panel.background=element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Probability density of the number of lineages through time")
  
  return (p)
}