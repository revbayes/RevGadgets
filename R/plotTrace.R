#' Plot trace
#'
#' Plots the posterior distributions of variables from trace file.
#'
#' Plots the posterior distributions of continuous variables from one or
#' multiple traces (as in, from multiple runs). Shaded regions under the curve
#' represent the 95\% credible interval.  If multiple traces are provided,
#' plotTrace() will plot each run independently as well as plot the combined
#' output. Note that for variables ith very different distributions, overlaying
#' the plots may result in illegible figures. In these cases, we recommend
#'  plotting each parameter separately.
#'
#'
#'
#' @param trace (list of data frames; no default) Name of a list of data frames,
#' such as produced by readTrace(). If the readTrace() output
#' contains multiple traces (such as from multiple runs), summarizeTrace() will
#' provide summaries for each trace individually, as well as the combined trace.
#'
#' @param color ("character"; "default") Colors for parameters. Defaults to
#' default RevGadgets colors. For non-default colors, provide a named vector of
#' length of the number of parameters.
#'
#' @param vars (character or character vector; NULL) The specific name(s) of
#' the variable(s) to be summarized.
#'
#' @param match (character; NULL) A string to match to a group of parameters.
#' For example, match = "er" will plot the variables "er[1]", "er[2]", "er[3]",
#' etc.. match will only work if your search string is followed by brackets in
#' one or more of the column names of the provided trace file. match = "er" will
#'only return the exchangeability parameters, but will not plot "Posterior".
#'
#' @return plotTrace() returns a list of the length of provided trace object,
#' plus one combined trace. Each element of the list contains a ggplot object
#' with plots of the provided parameters. These plots may be modified in
#' typical ggplot fashion.
#'
#' @examples
#'
#' \donttest{
#'
#' # example with quantitative parameters
#'
#' # download the example dataset to working directory
#' url_gtr <-
#'    "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.log"
#' dest_path_gtr <- "primates_cytb_GTR.log"
#' download.file(url_gtr, dest_path_gtr)
#'
#' # to run on your own data, change this to the path to your data file
#' file <- dest_path_gtr
#'
#' one_trace <- readTrace(paths = file)
#' plots <- plotTrace(trace = one_trace,
#'                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
#' plots[[1]]
#'
#' # add custom colors
#' plots <- plotTrace(trace = one_trace,
#'                    vars = c("pi[3]","pi[4]","pi[1]","pi[2]"),
#'                    color = c("pi[1]" = "green",
#'                              "pi[2]"= "red",
#'                              "pi[3]"= "blue",
#'                              "pi[4]"= "orange"))
#' plots[[1]]
#'
#' # make the same plot, using match
#' plots <- plotTrace(trace = one_trace, match = "pi")
#' plots[[1]]
#'
#' #' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_gtr)
#'
#' # plot some qualitative variables
#'
#' # download the example dataset to working directory
#' url_rj <- "https://revbayes.github.io/tutorials/intro/data/freeK_RJ.log"
#' dest_path_rj <- "freeK_RJ.log"
#' download.file(url_rj, dest_path_rj)
#'
#' file <- dest_path_rj
#' trace <- readTrace(path = file)
#'
#' plots <- plotTrace(trace = trace,
#'                    vars = c("prob_rate_12", "prob_rate_13",
#'                             "prob_rate_31", "prob_rate_32"))
#' plots[[1]]
#'
#' # with custom colors
#' plots <- plotTrace(trace = trace,
#'                    vars = c("prob_rate_12", "prob_rate_13",
#'                             "prob_rate_31", "prob_rate_32"),
#'                    color = c("prob_rate_12" = "green",
#'                              "prob_rate_13" = "red",
#'                              "prob_rate_31"= "blue",
#'                              "prob_rate_32" = "orange"))
#' plots[[1]]
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_rj)
#'
#' }
#'
#' @export


plotTrace <-
  function(trace,
           color = "default",
           vars = NULL,
           match = NULL) {
    # enforce argument matching
    if (is.list(trace) == FALSE)
      stop("trace should be a list of data frames")
    if (is.data.frame(trace[[1]]) == FALSE)
      stop("trace should be a list of data frames")
    if (is.character(color) == FALSE)
      stop ("color should be 'default' or valid color(s)")
    if (color[1] != "default" &
        any(.isColor(color) == FALSE))
      stop("node_color should be valid color(s)")
    if (color[1] == "default" & length(vars) > 12) {
      stop(
        paste0(
          length(vars),
          " states in dataset; please provide colors
          (default only can provide up to 12)"
        )
      )
    }
    if (color[1] != "default" & length(color) < length(vars)) {
      stop(
        paste0(
          "You provided fewer colors in node_color than states in your dataset.
          There are ",
          length(vars),
          " states and you provide ",
          length(color),
          " colors."
        )
      )
    }
    if (color[1] != "default" & length(color) > length(vars)) {
      stop(
        paste0(
          "You provided more colors in node_color than states in your dataset.
          There are ",
          length(vars),
          " states and you provide ",
          length(color),
          " colors."
        )
      )
    }
    if (is.null(vars) &
        is.null(match))
      stop("Either vars or match should be specified")
    if (!is.null(vars) &
        !is.null(match))
      stop("Only vars OR match should be specified")
    if (is.character(vars) == FALSE &
        !is.null(vars))
      stop("vars should be a character vector")
    if (is.character(match) == FALSE &
        !is.null(match))
      stop("match should be a character vector")

    # ensure variable names present in data frame
    if (!is.null(vars)) {
      if (any(vars %in% colnames(trace[[1]]) == FALSE) == TRUE) {
        stop(paste0(
          "The following variables you provided are not present in trace file:",
          paste0("\t", vars[!vars %in% colnames(trace[[1]])]),
          sep = "\n")
        )
      }
    }

    # find matching column names if using match
    if (!is.null(match)) {
      colnames <-
        unlist(lapply(strsplit(colnames(trace[[1]]), "\\["), function(x)
          x[1]))
      vars <- colnames(trace[[1]])[colnames %in% match]
      if (length(vars) == 0) {
        stop("match did not correspond to any column names in provided trace")
      }
    }


    plots <- list()
    #identify type of vars and split by quantitative vs. qualitative
    classes <- character()
    for (i in seq_len(length(vars))) {
      classes[i] <- class(trace[[1]][, vars[i]])
    }
    vars_quant <- vars[classes == "numeric"]
    vars_qual <- vars[classes != "numeric"]
    if (length(vars_quant) > 12) {
      stop("Please supply fewer than 12 quantitative variables")
    }
    if (length(vars_qual) > 12) {
      stop("Please supply fewer than 12 qualitative variables")
    }

    # set up colors
    if (color[1] == "default") {
      if (length(vars_quant) > 0) {
        col_vec_quant <- colFun(length(vars_quant))
        names(col_vec_quant) <- vars[classes == "numeric"]
      }
      if (length(vars_qual) > 0) {
        col_vec_qual <- colFun(length(vars_qual))
        names(col_vec_qual) <- vars[classes != "numeric"]
      }
    } else {
      col_vec_quant <- color[vars_quant]
      col_vec_qual <- color[vars_qual]
    }

    # make the quantitative plots
    if (length(vars_quant) > 0) {
      for (i in seq_len(length(trace))) {
        if (length(vars_quant) > 1) {
          # reshape data for plotting
          t <- trace[[i]][, vars_quant]
          t <- reshape::melt(t)
          colnames(t) <- c("Variable", "value")
          dfs <- list()
          for (k in seq_len(length(vars_quant))) {
            den <- density(trace[[i]][, vars_quant[k]])
            dfs[[k]] <-
              data.frame(Variable = vars_quant[k],
                         x = den$x,
                         y = den$y)
          }
          tt <- do.call("rbind", dfs)
          # calculate quantiles for all variables
          q_lows <- numeric()
          q_highs <- numeric()
          for (k in seq_len(length(vars_quant))) {
            q_lows[k] <-
              quantile(t[t$Variable == vars_quant[k], "value"], 0.025)
          }
          for (k in seq_len(length(vars_quant))) {
            q_highs[k] <-
              quantile(t[t$Variable == vars_quant[k], "value"], 0.975)
          }

          # plot densities, filling in the 95% credible interval
          plots[[i]] <- ggplot2::ggplot(t) +
            ggplot2::stat_density(
              ggplot2::aes(x = value, color = Variable),
              geom = "line",
              position = "identity"
            ) +
            ggplot2::xlab("Parameter value") +
            ggplot2::ylab("Density")
          for (k in seq_len(length(vars_quant))) {
            plots[[i]] <- plots[[i]] +
              ggplot2::geom_ribbon(
                data = subset(tt, Variable == vars_quant[k] &
                                x > q_lows[k] & x < q_highs[k]),
                ggplot2::aes(
                  x = x,
                  ymax = y,
                  ymin = 0
                ),
                fill = col_vec_quant[k],
                alpha = 0.5
              )
          }

          plots[[i]] <- plots[[i]] +
            ggplot2::scale_color_manual(values = col_vec_quant) +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = paste("Trace", i, sep = " ")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        } else if (length(vars_quant) == 1) {
          t <- data.frame(value = trace[[i]][, vars_quant])
          q_low <- quantile(t$value, 0.025)
          q_high <- quantile(t$value, 0.975)
          tt <-
            data.frame(x = density(t$value)$x, y = density(t$value)$y)

          # plot density, filling in the 95% credible interval
          plots[[i]] <- ggplot2::ggplot(t) +
            ggplot2::stat_density(
              ggplot2::aes(x = value),
              position = "identity",
              color = col_vec_quant[1],
              geom = "line"
            ) +
            ggplot2::xlab("Parameter value") +
            ggplot2::ylab("Density") +
            ggplot2::geom_ribbon(
              data = subset(tt, x > q_low & x < q_high),
              ggplot2::aes(
                x = x,
                ymax = y,
                ymin = 0
              ),
              fill = col_vec_quant[1]
            ) +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = paste0("Trace ", i, ": ", vars_quant)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }
      }
      # add combined trace plots if multiple traces provided
      if (length(trace) > 1) {
        if (length(vars_quant) > 1) {
          t <- do.call("rbind", trace)[, vars_quant]
          dfs <- list()
          for (k in seq_len(length(vars_quant))) {
            den <- density(t[, vars_quant[k]])
            dfs[[k]] <-
              data.frame(Variable = vars_quant[k],
                         x = den$x,
                         y = den$y)
          }
          tt <- do.call("rbind", dfs)
          t <- reshape::melt(t)
          colnames(t) <- c("Variable", "value")
          # calculate quantiles for all variables
          q_lows <- numeric()
          q_highs <- numeric()
          for (k in seq_len(length(vars_quant))) {
            q_lows[k] <-
              quantile(t[t$Variable == vars_quant[k], "value"], 0.025)
          }
          for (k in seq_len(length(vars_quant))) {
            q_highs[k] <-
              quantile(t[t$Variable == vars_quant[k], "value"], 0.975)
          }

          # plot densities, filling in the 95% credible interval
          plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
            ggplot2::stat_density(
              ggplot2::aes(x = value,
                           color = Variable),
              position = "identity",
              geom = "line"
            ) +
            ggplot2::xlab("Parameter value") +
            ggplot2::ylab("Density")
          for (k in seq_len(length(vars_quant))) {
            plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
              ggplot2::geom_ribbon(
                data = subset(tt, Variable == vars_quant[k] &
                                x > q_lows[k] & x < q_highs[k]),
                ggplot2::aes(
                  x = x,
                  ymax = y,
                  ymin = 0
                ),
                fill = col_vec_quant[k],
                alpha = 0.5
              )
          }

          plots[[length(trace) + 1]] <- plots[[length(trace) + 1]] +
            ggplot2::scale_color_manual(values = col_vec_quant) +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = "Combined Trace:")  +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        } else if (length(vars_quant) == 1) {
          t <- data.frame(value =  do.call("rbind", trace)[, vars_quant])
          q_low <- quantile(t$value, 0.025)
          q_high <- quantile(t$value, 0.975)
          tt <-
            data.frame(x = density(t$value)$x, y = density(t$value)$y)

          # plot density, filling in the 95% credible interval
          plots[[length(trace) + 1]] <- ggplot2::ggplot(t) +
            ggplot2::stat_density(
              ggplot2::aes(x = value),
              color = col_vec_quant[1],
              geom = "line",
              position = "identity"
            ) +
            ggplot2::geom_ribbon(
              data = subset(tt, x > q_low & x < q_high),
              ggplot2::aes(
                x = x,
                ymax = y,
                ymin = 0
              ),
              fill = col_vec_quant[1]
            ) +
            ggplot2::xlab("Parameter value") +
            ggplot2::ylab("Density") +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = paste0("Combined Trace: ", vars_quant)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }
      }
    }

    nplots <- length(plots)

    # make the qualitative plots
    if (length(vars_qual > 0)) {
      for (i in seq_len(length(trace))) {
        if (length(vars_qual) > 1) {
          # subset trace by qualitative variables
          t <- trace[[i]][, vars_qual]
          # calculate state probabilities and the credible set for all variables
          sp <- list()
          for (k in seq_len(length(vars_qual))) {
            # add in catch case for when table only 1 state visited
            table <-
              sort(table(t[, vars_qual[k]]) / nrow(t), decreasing = TRUE)
            if (is.table(table) == FALSE) {
              sp[[k]] <- data.frame(
                Var1 = names(table),
                Freq = 1,
                var = vars_qual[k],
                cred_set = TRUE
              )
            } else {
              sp_df <-
                data.frame(table,
                           var = vars_qual[k],
                           stringsAsFactors = FALSE)
              cs_table <-
                as.table(table[1:min(which((cumsum(table) >= 0.95) == TRUE))])
              cs_df <-
                data.frame(cs_table,
                           var = vars_qual[k],
                           stringsAsFactors = FALSE)
              sp[[k]] <-
                data.frame(sp_df, cred_set = sp_df$Freq %in% cs_df$Freq)
            }
          }
          state_probs <- do.call("rbind", sp)
          colnames(state_probs) <-
            c("State", "Probability", "Variable", "cred_set")
          state_probs$col <- state_probs$Variable
          state_probs$col[!state_probs$cred_set] <- "zzzzzzz"
          # reorder
          state_probs$Variable <- factor(state_probs$Variable,
                                         levels = unique(state_probs$Variable))
          state_probs$col <-
            factor(state_probs$col,
                   levels = c(levels(state_probs$Variable), "zzzzzzz"))
          state_probs$State <-
            as.character(levels(state_probs$State))[state_probs$State]
          # fill in 0 probabilities so all variables have data for all states
          state_probs <-
            tidyr::complete(
              data = state_probs,
              State,
              Variable,
              fill = list(
                probability = 0,
                cred_set = FALSE,
                col = "zzzzzzz"
              )
            )

          # plot the variables
          plots[[nplots + i]] <-
            ggplot2::ggplot(
              state_probs,
              ggplot2::aes(
                color = Variable,
                fill = col,
                y = Probability,
                x = State
              )
            ) +
            ggplot2::geom_bar(
              position = ggplot2::position_dodge2(preserve = "single"),
              stat = "identity"
            ) +
            ggplot2::theme_bw() +
            ggplot2::scale_color_manual(values = col_vec_qual) +
            ggplot2::scale_fill_manual(values = c(col_vec_qual, "#ffffff"),
                                       guide = "none") +
            ggplot2::guides(color =
                              ggplot2::guide_legend(override.aes =
                                                      list(fill =
                                                             col_vec_qual))) +
            ggplot2::ggtitle(label = paste0("Trace ", i)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }

        if (length(vars_qual) == 1) {
          # subset trace by qual variable
          t <- trace[[i]][, vars_qual]

          # calculate credible set and rearrange data for plotting
          state_probs <- sort(table(t) / length(t), decreasing = TRUE)
          data_full <- data.frame(state_probs)
          colnames(data_full) <- c("State", "Probability")
          data_full$State <-
            as.character(levels(data_full$State))[data_full$State]
          data_sig <-
            data_full[1:min(which((
              cumsum(data_full$Probability) >= 0.95
            ) == TRUE)), ]

          #plot
          plots[[nplots + i]] <-
            ggplot2::ggplot(data_full,
                            ggplot2::aes(x = State, y = Probability)) +
            ggplot2::geom_bar(stat = "identity",
                              color = col_vec_qual[1],
                              fill = "white") +
            ggplot2::geom_bar(
              data = data_sig,
              stat = "identity",
              ggplot2::aes(x = State, y = Probability),
              color = col_vec_qual[1],
              fill = col_vec_qual[1]
            ) +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = paste0("Trace ", i, ": ", vars_qual)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }
      }
      # add combined trace for multiple traces
      nplots <- length(plots)  #reset nplots
      if (length(trace) > 1) {
        if (length(vars_qual) > 1) {
          # subset trace by qualitative variables
          t <- do.call("rbind", trace)[, vars_qual]
          # calculate state probabilities and the credible set for all variables
          sp <- list()
          for (k in seq_len(length(vars_qual))) {
            table <- sort(table(t[, vars_qual[k]]) / nrow(t), decreasing = TRUE)
            if (is.table(table) == FALSE) {
              sp[[k]] <- data.frame(
                Var1 = names(table),
                Freq = 1,
                var = vars_qual[k],
                cred_set = TRUE
              )
            } else {
              sp_df <-
                data.frame(table,
                           var = vars_qual[k],
                           stringsAsFactors = FALSE)
              cs_table <-
                as.table(table[1:min(which((cumsum(table) >= 0.95) == TRUE))])
              cs_df <-
                data.frame(cs_table,
                           var = vars_qual[k],
                           stringsAsFactors = FALSE)
              sp[[k]] <-
                data.frame(sp_df, cred_set = sp_df$Freq %in% cs_df$Freq)
            }
          }
          state_probs <- do.call("rbind", sp)
          colnames(state_probs) <-
            c("State", "Probability", "Variable", "cred_set")
          state_probs$col <- state_probs$Variable
          state_probs$col[!state_probs$cred_set] <- "zzzzzzz"
          state_probs$State <-
            as.character(levels(state_probs$State))[state_probs$State]
          # fill in 0 probabilities so all variables have data for all states
          state_probs <-
            tidyr::complete(
              data = state_probs,
              State,
              Variable,
              fill = list(
                Probability = 0,
                cred_set = FALSE,
                col = "zzzzzzz"
              )
            )
          # plot the variables
          plots[[nplots + 1]] <-
            ggplot2::ggplot(
              state_probs,
              ggplot2::aes(
                color = Variable,
                fill = col,
                y = Probability,
                x = State
              )
            ) +
            ggplot2::geom_bar(
              position = ggplot2::position_dodge2(preserve = "single"),
              stat = "identity"
            ) +
            ggplot2::theme_bw() +
            ggplot2::scale_color_manual(values = col_vec_qual) +
            ggplot2::scale_fill_manual(values = c(col_vec_qual, "#ffffff"),
                                       guide = "none") +
            ggplot2::guides(color =
                              ggplot2::guide_legend(override.aes =
                                                       list(fill =
                                                              col_vec_qual))) +
            ggplot2::ggtitle(label = paste0("Combined Trace")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        } else if (length(vars_qual) == 1) {
          #combine all traces
          #subset trace by qual variable
          t <- do.call("rbind", trace)[, vars_qual]

          #calculate credible set and rearrange data for plotting
          state_probs <- sort(table(t) / length(t), decreasing = TRUE)
          data_full <- data.frame(state_probs)
          colnames(data_full) <- c("State", "Probability")
          data_full$State <-
            as.character(levels(data_full$State))[data_full$State]
          data_sig <-
            data_full[1:min(which((
              cumsum(data_full$Probability) >= 0.95
            ) == TRUE)), ]

          #plot
          plots[[nplots + 1]] <-
            ggplot2::ggplot(data_full,
                            ggplot2::aes(x = State, y = Probability)) +
            ggplot2::geom_bar(stat = "identity",
                              color = col_vec_qual[1],
                              fill = "white") +
            ggplot2::geom_bar(
              data = data_sig,
              stat = "identity",
              ggplot2::aes(x = State, y = Probability),
              color = col_vec_qual[1],
              fill = col_vec_qual[1]
            ) +
            ggplot2::theme_bw() +
            ggplot2::ggtitle(label = paste0("Combined Trace: ", vars_qual)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        }
      }
    }

    return(plots)
  }
