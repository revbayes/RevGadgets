#' Plot Validation Results from RevBayes Analyses
#'
#' @description
#' Creates a ggplot visualization of RevBayes validation results. It plots the observed coverage
#' probability against the expected probability (the HPD width). A
#' perfectly calibrated model should yield a diagonal line from (0,0) to (1,1).
#'
#' @param validation_results (list; no default) A named list of data frames as returned by the
#'  `processValidation` function, where each data frame corresponds to a parameter and contains
#'  the HPD width, the count of simulations where the true value was within the HPD interval,
#'  the total count and the coverage frequency.
#'
#' @param parameter (character string; optional) The name of a single parameter to plot from the
#'   'validation_results' list. If 'NULL' (the default), the function will generate a list of plots for all
#'   parameters found in the input.
#'
#' @return A ggplot object if a single parameter is specified, or a named list of ggplot objects
#'   if multiple parameters are present in the input list.
#'
#' @examples
#' \dontrun{
#' # First, generate the validation results
#' validation_output <- processValidation(path = "output_my_analysis", n_reps = 100)
#'
#' # Plot the results for a single parameter
#' p_alpha <- plotValidation(validation_output, parameter = "alpha")
#'
#' # Display the plot
#' p_alpha
#'
#' # Customize and save it
#' p_alpha <- p_alpha + ggplot2::labs(title = "Validation for Alpha Parameter")
#' ggplot2::ggsave("alpha_validation.png", p_alpha, width = 5, height = 5)
#' }
#'
#' @export
#' @import ggplot2
#'
plotValidation <- function(validation_results, parameter = NULL) {
  # input validation
  if (!is.list(validation_results) || length(validation_results) == 0) {
    stop("Argument 'validation_results' must be a non-empty list, ideally as returned by the processValidation() function.")
  }

  # determine which parameters to plot
  if (is.null(parameter)) {
    # plot all the parameters found in the list
    parameters_to_plot <- names(validation_results)
  } else {
    # plot only the specified parameter
    if (!parameter %in% names(validation_results)) {
      stop(paste("Parameter", parameter, "not found in 'validation_results'. Available parameters are:", paste(names(validation_results), collapse = ", ")))
    }
    parameters_to_plot <- parameter
  }

  plot_list <- list()

  # loop through the selected parameters
  for (param in parameters_to_plot) {
    coverage_data <- validation_results[[param]]

    # check if the data frame has the required columns
    required_cols <- c("hpd_width", "freq")
    if (!all(required_cols %in% colnames(coverage_data))) {
      warning("Skipping parameter '", param, "' due to missing columns in its data frame.")
      next
    }

    # the plot
    p <- ggplot(coverage_data, aes(x = .data$hpd_width, y = .data$freq)) +
      geom_bar(stat = "identity", colour = "lightgray", fill = "lightgray") +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
        linetype = "dashed", linewidth = 0.8, color = "black"
      ) +
      labs(
        x = "HPD Credible Interval Width",
        y = "Coverage Probability",
        title = param
      ) +
      theme_classic(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    plot_list[[param]] <- p
  }

  # return the plot object(s)
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}
