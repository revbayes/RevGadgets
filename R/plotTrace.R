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
    #identify type of vars and split by quantitative vs. qualitative
    classes <- character()
    for (i in 1:length(vars)) {
      classes[i] <- class(trace[[1]][,vars[i]])
    }
    vars_quant <- vars[classes == "numeric"]
    vars_qual <- vars[classes != "numeric"]
    if (length(vars_quant) > 12) {stop("Please supply fewer than 12 quantitative variables")}
    if (length(vars_qual) > 12) {stop("Please supply fewer than 12 qualitative variables")}

    # make the quantitative plots
    if (length(vars_quant > 0)) {
    for(i in 1:length(trace)){
      if (length(vars_quant) > 1) {
        # reshape data for plotting
        t <- trace[[i]][,vars_quant]
        t <- reshape::melt(t)
        dfs <- list()
        for (k in 1:length(vars_quant)) {
          den<- density(trace[[i]][,vars_quant[k]])
          dfs[[k]] <- data.frame(variable=vars_quant[k], x = den$x, y = den$y)
        }
        tt <- do.call("rbind", dfs)
        # calculate quantiles for all variables
        q_lows <- numeric()
        q_highs <- numeric()
        for (k in 1:length(vars_quant)) {q_lows[k] <- quantile(t[t$variable == vars_quant[k],"value"],0.025)}
        for (k in 1:length(vars_quant)) {q_highs[k] <- quantile(t[t$variable == vars_quant[k],"value"],0.975)}

        # plot densities, filling in the 95% credible interval
        plots[[i]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value, color = variable),
                                geom="line", position = "identity")
        for (k in 1:length(vars_quant)) {
          plots[[i]] <- plots[[i]] +
            ggplot2::geom_ribbon(data = subset(tt, variable == vars_quant[k] & x > q_lows[k] & x < q_highs[k]),
                                        ggplot2::aes(x=x,ymax=y,ymin=0),
                                        fill = .colFun(length(vars_quant))[k], alpha=0.5)
        }

        plots[[i]] <- plots[[i]] +
          ggplot2::scale_color_manual(values = .colFun(length(vars_quant))) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste("Trace",i, sep = " ") ) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[i]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value,
        #                                                     fill = variable,
        #                                                     color = variable)) +
        #              ggplot2::geom_density(alpha = 0.5) +
        #              ggplot2::scale_fill_manual(values = .colFun(length(vars_quant))) +
        #              ggplot2::scale_color_manual(values = .colFun(length(vars_quant))) +
        #              ggthemes::theme_few() +
        #              ggplot2::ggtitle(label = paste("Trace",i, sep = " ") ) +
        #              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


      } else if (length(vars_quant) == 1) {
        t <- data.frame(value = trace[[i]][,vars_quant])
        q_low <- quantile(t$value, 0.025)
        q_high <- quantile(t$value, 0.975)
        tt <- data.frame(x = density(t$value)$x, y = density(t$value)$y)

        # plot density, filling in the 95% credible interval
        plots[[i]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value),
                                position = "identity",
                                color = .colFun(1),geom="line") +
          ggplot2::geom_ribbon(data = subset(tt, x > q_low & x < q_high),
                               ggplot2::aes(x=x,ymax=y,ymin=0),
                               fill = .colFun(1)) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste0("Trace ",i,": ",vars_quant)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[i]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value)) +
        #  ggplot2::geom_density(fill = .colFun(1), color = .colFun(1)) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = paste0("Trace ",i,": ",vars_quant)) +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      }
    }
    # add combined trace plots if multiple traces provided
    if (length(trace) > 1){
      if (length(vars_quant) > 1) {
        t <- do.call("rbind", trace)[,vars_quant]
        dfs <- list()
        for (k in 1:length(vars_quant)) {
          den<- density(t[,vars_quant[k]])
          dfs[[k]] <- data.frame(variable=vars_quant[k], x = den$x, y = den$y)
        }
        tt <- do.call("rbind", dfs)
        t <- reshape::melt(t)

        # calculate quantiles for all variables
        q_lows <- numeric()
        q_highs <- numeric()
        for (k in 1:length(vars_quant)) {q_lows[k] <- quantile(t[t$variable == vars_quant[k],"value"],0.025)}
        for (k in 1:length(vars_quant)) {q_highs[k] <- quantile(t[t$variable == vars_quant[k],"value"],0.975)}

        # plot densities, filling in the 95% credible interval
        plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value,
                                             color = variable),
                                position = "identity",
                                geom="line")
        for (k in 1:length(vars_quant)) {
          plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
            ggplot2::geom_ribbon(data = subset(tt, variable == vars_quant[k] & x > q_lows[k] & x < q_highs[k]),
                                 ggplot2::aes(x=x,ymax=y,ymin=0),
                                 fill = .colFun(length(vars_quant))[k], alpha=0.5)
        }

        plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
          ggplot2::scale_color_manual(values = .colFun(length(vars_quant))) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = "Combined Trace:")  +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        # old version with no differential fill for credible interval
        #plots[[length(trace) + 1]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value,
        #                                                                     fill = variable,
        #                                                                     color = variable)) +
        #  ggplot2::geom_density(alpha = 0.5) +
        #  ggplot2::scale_fill_manual(values = .colFun(length(vars_quant))) +
        #  ggplot2::scale_color_manual(values = .colFun(length(vars_quant))) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = "Combined Trace:") +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      } else if (length(vars_quant) == 1) {
        t <- data.frame(value =  do.call("rbind", trace)[,vars_quant])
        #t <- data.frame(value = trace[[i]][,vars_quant])
        q_low <- quantile(t$value, 0.025)
        q_high <- quantile(t$value, 0.975)
        tt <- data.frame(x = density(t$value)$x, y = density(t$value)$y)

        # plot density, filling in the 95% credible interval
        plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
          ggplot2::stat_density(ggplot2::aes(x = value),
                                color = .colFun(1),geom="line",
                                position = "identity") +
          ggplot2::geom_ribbon(data = subset(tt, x > q_low & x < q_high),
                               ggplot2::aes(x=x,ymax=y,ymin=0),
                               fill = .colFun(1)) +
          ggthemes::theme_few() +
          ggplot2::ggtitle(label = paste("Combined Trace")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


        # old version with no differential fill for credible interval
        #plots[[length(trace) + 1]] <- ggplot2::ggplot(data = t, ggplot2::aes(x = value)) +
        #  ggplot2::geom_density(fill = .colFun(1), color = .colFun(1)) +
        #  ggthemes::theme_few() +
        #  ggplot2::ggtitle(label = paste("Combined Trace:", vars_quant)) +
        #  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      }
    }
    }

     nplots <- length(plots)

    if (length(vars_qual > 0)) {
      for(i in 1:length(trace)){
        if (length(vars_qual) > 1) {
          #subset trace by qualitative variables
          t <- trace[[i]][,vars_qual]
          #calculate state probabilities and the credible set for all variables
          sp <- list()
          for (k in 1:length(vars_qual)) {
            table <- sort(table(t[,vars_qual[k]])/nrow(t), decreasing = TRUE)
            sp_df <- data.frame(table, var = vars_qual[k], stringsAsFactors = F)
            cs_table <- as.table(table[1:min(which((cumsum(table) >= 0.95) == TRUE))])
            cs_df <- data.frame(cs_table,var = vars_qual[k], stringsAsFactors = F)
            sp[[k]] <- data.frame(sp_df, cred_set = sp_df$Freq %in% cs_df$Freq) # will this always work?
          }
          state_probs <- do.call("rbind", sp)
          colnames(state_probs) <- c("state","probability","variable", "cred_set")
          state_probs$col <- state_probs$variable
          state_probs$col[!state_probs$cred_set] <- "zzzzzzz"  #cheat so it's always last and will match up with the white hex code
          state_probs$state <- as.integer(state_probs$state)
          #reorder
          state_probs$variable <- factor(state_probs$variable,
                                         levels = unique(state_probs$variable))
          state_probs$col <- factor(state_probs$col,
                                         levels = c(levels(state_probs$variable), "zzzzzzz"))

          #fill in 0 probabilities such that all variables have data for all states
          state_probs <- tidyr::complete(data = state_probs, state, variable,
                                  fill = list(probability = 0,
                                              cred_set = FALSE,
                                              col = "zzzzzzz"))

          # plot the variables
          plots[[nplots + i]] <-
            ggplot2::ggplot(state_probs,
                            ggplot2::aes(color = variable, fill = col,
                                         y = probability, x = state)) +
            ggplot2::geom_bar(position = position_dodge2(preserve = "single"),
                              stat = "identity") +
            ggthemes::theme_few() +
            ggplot2::scale_color_manual(values = .colFun(length(vars_qual))) +
            ggplot2::scale_fill_manual(values = c(.colFun(length(vars_qual)), "#ffffff"), guide = FALSE) +
            ggplot2::guides(color = guide_legend(override.aes=list(fill=.colFun(length(vars_qual))))) +
            ggplot2::ggtitle(label = paste0("Trace ", i)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }

        if (length(vars_qual) == 1) {
          #subset trace by qual variable
          t <- trace[[i]][,vars_qual]

          #calculate credible set and rearrange data for plotting
          state_probs <- sort(table(t)/length(t), decreasing = TRUE)
          data_sig <- data.frame(table(state_probs[1:min(which((cumsum(state_probs) >= 0.95) == TRUE))]))
          data_full <- data.frame(state_probs)
          colnames(data_full) <- c("state", "probability")
          colnames(data_sig) <- c("probability", "state")
          data_sig$probability <- as.numeric(levels(data_sig$probability)[data_sig$probability])

          #plot
          plots[[nplots + i]] <- ggplot2::ggplot(data_full, ggplot2::aes(x = state, y = probability)) +
                                ggplot2::geom_bar(stat = "identity", color = .colFun(1), fill = "white") +
                                ggplot2::geom_bar(data = data_sig, stat = "identity",
                                                  ggplot2::aes(x = state, y = probability),
                                        color = .colFun(1), fill = .colFun(1)) +
                                ggthemes::theme_few() +
                                ggplot2::ggtitle(label = paste0("Trace ",i,": ",vars_qual)) +
                                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))






        }
      }
      # add combined trace for multiple traces
      nplots <- length(plots)  #reset nplots
      if (length(trace) > 1){
        if (length(vars_qual) > 1) {
          #subset trace by qualitative variables
          t <- do.call("rbind", trace)[,vars_qual]
          #calculate state probabilities and the credible set for all variables
          sp <- list()
          for (k in 1:length(vars_qual)) {
            table <- sort(table(t[,vars_qual[k]])/nrow(t), decreasing = TRUE)
            sp_df <- data.frame(table, var = vars_qual[k], stringsAsFactors = F)
            cs_table <- as.table(table[1:min(which((cumsum(table) >= 0.95) == TRUE))])
            cs_df <- data.frame(cs_table,var = vars_qual[k], stringsAsFactors = F)
            sp[[k]] <- data.frame(sp_df, cred_set = sp_df$Freq %in% cs_df$Freq) # will this always work?
          }
          state_probs <- do.call("rbind", sp)
          colnames(state_probs) <- c("state","probability","variable", "cred_set")
          state_probs$col <- state_probs$variable
          state_probs$col[!state_probs$cred_set] <- "zzzzzzz"  #cheat so it's always last and will match up with the white hex code

          #fill in 0 probabilities such that all variables have data for all states
          state_probs <- tidyr::complete(data = state_probs, state, variable,
                                         fill = list(probability = 0,
                                                     cred_set = FALSE,
                                                     col = "zzzzzzz"))
          # plot the variables
          plots[[nplots + 1]] <-
            ggplot2::ggplot(state_probs,
                            ggplot2::aes(color = variable, fill = col,
                                         y = probability, x = state)) +
            ggplot2::geom_bar(position = position_dodge2(preserve = "single"),
                              stat = "identity") +
            ggthemes::theme_few() +
            ggplot2::scale_color_manual(values = .colFun(length(vars_qual))) +
            ggplot2::scale_fill_manual(values = c(.colFun(length(vars_qual)), "#ffffff"), guide = FALSE) +
            ggplot2::guides(color = ggplot2::guide_legend(override.aes=list(fill=.colFun(length(vars_qual))))) +
            ggplot2::ggtitle(label = paste0("Combined Trace")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        } else if (length(vars_qual) == 1) {

          #combine all traces
          #subset trace by qual variable
          t <- do.call("rbind", trace)[,vars_qual]

          #calculate credible set and rearrange data for plotting
          state_probs <- sort(table(t)/length(t), decreasing = TRUE)
          data_sig <- data.frame(table(state_probs[1:min(which((cumsum(state_probs) >= 0.95) == TRUE))]))
          data_full <- data.frame(state_probs)
          colnames(data_full) <- c("state", "probability")
          colnames(data_sig) <- c("probability", "state")
          data_sig$probability <- as.numeric(levels(data_sig$probability)[data_sig$probability])

          #plot
          plots[[nplots + 1]] <- ggplot2::ggplot(data_full, ggplot2::aes(x = state, y = probability)) +
            ggplot2::geom_bar(stat = "identity", color = .colFun(1), fill = "white") +
            ggplot2::geom_bar(data = data_sig, stat = "identity",
                              aes(x = state, y = probability),
                              color = .colFun(1), fill = .colFun(1)) +
            ggthemes::theme_few() +
            ggplot2::ggtitle(label = paste0("Combined Trace: ",vars_qual)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        }
        }
    }

  return(plots)
}