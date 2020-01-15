#' Plot trace
#'
#' Plots the posterior distributions of variables from trace file.
#'
#' Plots the posterior distributions of continuous variables from one or
#' multiple traces (as in, from multiple runs). Shaded regions under the curve
#' represent the 95\% credible interval.  If multiple traces are provided,
#' plotTrace() will plot each run independently as well as plot the combined output.
#' Note that for variables ith very different distributions, overlaying the
#' plots may result in illegible figures. In these cases, we recommend plotting
#' each parameter separately.
#'
#'
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), summarizeTrace() will provide
#' summaries for each trace individually, as well as the combined trace.
#'
#' @param vars (character or character vector; NULL) The specific name(s) of the variable(s)
#' to be summarized.
#'
#' @param match (character; NULL) A string to match to a group of parameters. For example,
#' match = "er" will plot the variables "er[1]", "er[2]", "er[3]", etc.. match
#' will only work if your search string is followed by brackets in one or more of
#'  the column names of the provided trace file. match = "er" will only return the
#'  exchangeability parameters, but will not plot "Posterior".
#'
#' @return plotTrace() returns a list of the length of provided trace object, plus one
#' combined trace. Each element of the list contains a ggplot object with plots of
#' the provided parameters. These plots may be modified in typical ggplot fashion.
#'
#' @examples
#'
#' \dontrun{
#'
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_covariotide.p", package="RevGadgets")
#' one_trace <- readTrace(paths = file)
#' plots <- plotTrace(trace = one_trace,
#'                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
#' plots[[1]]
#' # add custom colors
#' plots[[1]] + scale_fill_manual(values = c("green","red","blue","orange"))
#'
#' # make the same plot, using match
#'
#' plots <- plotTrace(trace = one_trace, match = "pi")
#'
#' }
#'
#' @export


plotTrace <- function(trace, vars = NULL, match = NULL) {

  # enforce argument matching
  if (is.list(trace) == FALSE) stop("trace should be a list of data frames")
  if (is.data.frame(trace[[1]]) == FALSE) stop("trace should be a list of data frames")
  if (is.null(vars) & is.null(match)) stop("Either vars or match should be specified")
  if (!is.null(vars) & !is.null(match)) stop("Only vars OR match should be specified")
  if (is.character(vars) == FALSE & !is.null(vars)) stop("vars should be a character vector")
  if (is.character(match) == FALSE & !is.null(match)) stop("match should be a character vector")

  # ensure variable names present in data frame
  if (!is.null(vars)){
    if (any(vars %in% colnames(trace[[1]]) == FALSE) == TRUE) {
      cat("The following variables you provided are not present in trace file:",
          paste0("\t", vars[!vars %in% colnames(trace[[1]])]), sep = "\n")
      stop("oops!")
    }
  }
  # find matching column names if using match
  if (!is.null(match)) {
    vars <- colnames(trace[[1]])[grep(paste0("(",match,")(\\[)"),colnames(trace[[1]]))]
    if (length(vars) == 0) {stop("match did not correspond to any column names in provided trace")}
  }

    plots <- list()
    for(i in 1:length(trace)){
      if (length(vars) > 1) {

        # reshape data for plotting
        t <- trace[[i]][,vars]
        t <- reshape::melt(t)
        dfs <- list()
        for (k in 1:length(vars)) {
          den<- density(trace[[i]][,vars[k]])
          dfs[[k]] <- data.frame(variable=vars[k], x = den$x, y = den$y)
        }
        tt <- do.call("rbind", dfs)
        # calculate quantiles for all variables
        q_lows <- numeric()
        q_highs <- numeric()
        for (k in 1:length(vars)) {q_lows[k] <- quantile(t[t$variable == vars[k],"value"],0.025)}
        for (k in 1:length(vars)) {q_highs[k] <- quantile(t[t$variable == vars[k],"value"],0.975)}

        # plot densities, filling in the 95% credible interval
        plots[[i]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value,
                                             color = variable),
                                geom="line", position="dodge")
        for (k in 1:length(vars)) {
          plots[[i]] <- plots[[i]] +
            ggplot2::geom_ribbon(data = subset(tt, variable == vars[k] & x > q_lows[k] & x < q_highs[k]),
                                        ggplot2::aes(x=x,ymax=y,ymin=0),
                                        fill = colFun(length(vars))[k], alpha=0.5)
        }

        plots[[i]] <- plots[[i]] +
          ggplot2::scale_color_manual(values = colFun(length(vars))) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste("Trace",i, sep = " ") ) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[i]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value,
        #                                                     fill = variable,
        #                                                     color = variable)) +
        #              ggplot2::geom_density(alpha = 0.5) +
        #              ggplot2::scale_fill_manual(values = colFun(length(vars))) +
        #              ggplot2::scale_color_manual(values = colFun(length(vars))) +
        #              ggthemes::theme_few() +
        #              ggplot2::ggtitle(label = paste("Trace",i, sep = " ") ) +
        #              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


      } else if (length(vars) == 1) {
        t <- data.frame(value = trace[[i]][,vars])
        q_low <- quantile(t$value, 0.025)
        q_high <- quantile(t$value, 0.975)
        tt <- data.frame(x = density(t$value)$x, y = density(t$value)$y)

        # plot density, filling in the 95% credible interval
        plots[[i]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value),
                                color = colFun(1),geom="line", position="dodge") +
          ggplot2::geom_ribbon(data = subset(tt, x > q_low & x < q_high),
                               ggplot2::aes(x=x,ymax=y,ymin=0),
                               fill = colFun(1)) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste0("Trace ",i,": ",vars)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[i]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value)) +
        #  ggplot2::geom_density(fill = colFun(1), color = colFun(1)) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = paste0("Trace ",i,": ",vars)) +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      }
    }
    # add combined trace plots if multiple traces provided
    if (length(trace) > 1){
      if (length(vars) > 1) {
        t <- do.call("rbind", trace)[,vars]
        t <- reshape::melt(t)
        dfs <- list()
        for (k in 1:length(vars)) {
          den<- density(trace[[i]][,vars[k]])
          dfs[[k]] <- data.frame(variable=vars[k], x = den$x, y = den$y)
        }
        tt <- do.call("rbind", dfs)
        # calculate quantiles for all variables
        q_lows <- numeric()
        q_highs <- numeric()
        for (k in 1:length(vars)) {q_lows[k] <- quantile(t[t$variable == vars[k],"value"],0.025)}
        for (k in 1:length(vars)) {q_highs[k] <- quantile(t[t$variable == vars[k],"value"],0.975)}

        # plot densities, filling in the 95% credible interval
        plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value,
                                             color = variable),
                                geom="line", position="dodge")
        for (k in 1:length(vars)) {
          plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
            ggplot2::geom_ribbon(data = subset(tt, variable == vars[k] & x > q_lows[k] & x < q_highs[k]),
                                 ggplot2::aes(x=x,ymax=y,ymin=0),
                                 fill = colFun(length(vars))[k], alpha=0.5)
        }

        plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
          ggplot2::scale_color_manual(values = colFun(length(vars))) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = "Combined Trace:")  +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[length(trace) + 1]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value,
        #                                                                     fill = variable,
        #                                                                     color = variable)) +
        #  ggplot2::geom_density(alpha = 0.5) +
        #  ggplot2::scale_fill_manual(values = colFun(length(vars))) +
        #  ggplot2::scale_color_manual(values = colFun(length(vars))) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = "Combined Trace:") +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      } else if (length(vars) == 1) {
        t <- data.frame(value =  do.call("rbind", trace)[,vars])
        #t <- data.frame(value = trace[[i]][,vars])
        q_low <- quantile(t$value, 0.025)
        q_high <- quantile(t$value, 0.975)
        tt <- data.frame(x = density(t$value)$x, y = density(t$value)$y)

        # plot density, filling in the 95% credible interval
        plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value),
                                color = colFun(1),geom="line", position="dodge") +
          ggplot2::geom_ribbon(data = subset(tt, x > q_low & x < q_high),
                               ggplot2::aes(x=x,ymax=y,ymin=0),
                               fill = colFun(1)) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste("Combined Trace:")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


        # old version with no differential fill for credible interval
        #plots[[length(trace) + 1]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value)) +
        #  ggplot2::geom_density(fill = colFun(1), color = colFun(1)) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = paste("Combined Trace:", vars)) +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      }
    }

  return(plots)
}
