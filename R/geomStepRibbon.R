#' plot geom stepribbon for diversification rates
#'
#' Modified from \code{RmcdrPlugin.KMggplot2} step ribbon plots.
#'
#' \code{geom_stepribbon} is an extension of the \code{geom_ribbon}, and
#' is optimized for Kaplan-Meier plots with pointwise confidence intervals
#' or a confidence band.
#'
#' @seealso
#'   \code{\link[ggplot2]{geom_ribbon}} \code{geom_stepribbon}
#'   inherits from \code{geom_ribbon}.
#'   \code{geom_stepribbon} is modified from
#'   \code{RcmdrPlugin.KMggplot2::geom_stepribbon}.
#'
#' @inheritParams ggplot2::geom_ribbon
#'
#' @rdname geom_stepribbon
#' @examples
#'
#' huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
#'
#' h <- ggplot2::ggplot(huron, ggplot2::aes(year))
#'
#' h + geom_stepribbon(ggplot2::aes(ymin = level - 1, ymax = level + 1),
#'                     fill = "grey70") +
#'     ggplot2::geom_step(ggplot2::aes(y = level))
#'
#' # contrast ggplot2::geom_ribbon with geom_stepribbon:
#' h + ggplot2::geom_ribbon(ggplot2::aes(ymin = level - 1, ymax = level + 1),
#'                          fill = "grey70") +
#'     ggplot2::geom_line(ggplot2::aes(y = level))
#'
#' @importFrom ggplot2 layer GeomRibbon
#'
#' @export
#'
geom_stepribbon <- function(mapping     = NULL,
                            data        = NULL,
                            stat        = "identity",
                            position    = "identity",
                            na.rm       = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE,
                            ...) {
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ...)
  )
}

#' @rdname geom_stepribbon
#' @format NULL
#' @usage NULL
#' @export
GeomStepribbon <- ggplot2::ggproto(
  "GeomStepribbon",
  GeomRibbon,

  extra_params = c("na.rm"),

  draw_group = function(data,
                        panel_scales,
                        coord,
                        direction = "vh",
                        include_final = FALSE,
                        na.rm = FALSE) {
    n <- nrow(data)
    data <- as.data.frame(data)[order(data$x),]

    if (direction == "vh") {
      xs <- rep(1:n, each = 2)[-2 * n]
      ys <- c(1, rep(2:n, each = 2))
    } else if (direction == "hv") {
      xs <- c(1, rep(2:n, each = 2))
      ys <- rep(1:n, each = 2)[-2 * n]
    } else {
      abort("Parameter `direction` is invalid.")
    }
    if (!include_final) {
      xs <- tail(xs, n = -1)
      ys <- tail(ys, n = -1)
    }
    x <- data$x[xs]
    ymin <- data$ymin[ys]
    ymax <- data$ymax[ys]
    # recover()
    data_attr <-
      data[xs, setdiff(names(data), c("x", "ymin", "ymax"))]
    #data <-
    #  .new_data_frame(c(list(
    #    x = x, ymin = ymin, ymax = ymax
    #  ), data_attr))
    data <- cbind(data.frame(x = x, ymin = ymin, ymax = ymax),
                   data_attr)

    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)

  }

)
