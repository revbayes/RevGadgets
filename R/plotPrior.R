
#' @export
plotPrior <- function(distribution, col=colFun(1)) {

  # get the distribution
  distn <- .getPrior(distribution)

  # plot the distribution
  pp <- ggplot2::ggplot(data.frame()) +
    ggplot2::stat_function(fun = distn$fun, col = col) +
    ggplot2::stat_function(fun = distn$fun, xlim = distn$hpd, geom="area", fill=col, alpha=0.5) +
    ggplot2::xlim(distn$min, distn$max) +
    ggplot2::theme_few() +
    ggplot2::labs(x = "x", y = "prior probability of x", title = distribution) +


  # return the plot
  return(pp)

}



#' @export
.getPrior <- function(distribution) {

  # parse the distribution string
  dist_string <- parse(text = distribution)[[1]]

  # evaluate it
  dist <- eval(dist_string)

  # return the distribution object
  return(dist)

}

dnExponential <- function(rate = 1, offset = 0) {
  min  <- offset + 0
  max  <- offset + qexp(0.999, rate = rate)
  hpd  <- offset + qexp(c(0.025, 0.975), rate = rate)
  fun  <- function(x) dexp(x - offset, rate = rate)
  dist <- list(min = min,
               max = max,
               hpd = hpd,
               fun = fun)
  return(dist)
}
dnExp <- dnExponential

dnGamma <- function(shape = 1, rate = 1) {
  min  <- 0
  max  <- qgamma(0.999, shape = shape, rate = rate)
  hpd  <- qgamma(c(0.025, 0.975), shape = shape, rate = rate)
  fun  <- function(x) dgamma(x, shape = shape, rate = rate)
  dist <- list(min = min,
               max = max,
               hpd = hpd,
               fun = fun)
  return(dist)
}











#' @export
plotPrior <- function(distribution, col=colFun(1)) {

  # get the distribution
  distn <- .getPrior(distribution)

  # plot the distribution
  pp <- ggplot2::ggplot(data.frame()) +
    ggplot2::stat_function(fun = distn$fun, col = col) +
    ggplot2::stat_function(fun = distn$fun, xlim = distn$hpd, geom="area", fill=col, alpha=0.5) +
    ggplot2::xlim(distn$min, distn$max) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "x", y = "prior probability of x", title = distribution) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))

  # return the plot
  return(pp)

}



#' @export
.getPrior <- function(distribution) {

  # parse the distribution string
  dist_string <- parse(text = distribution)[[1]]

  # evaluate it
  dist <- eval(dist_string)

  # return the distribution object
  return(dist)

}

dnExponential <- function(rate = 1, offset = 0) {
  min  <- offset + 0
  max  <- offset + qexp(0.999, rate = rate)
  hpd  <- offset + qexp(c(0.025, 0.975), rate = rate)
  fun  <- function(x) dexp(x - offset, rate = rate)
  dist <- list(min = min,
               max = max,
               hpd = hpd,
               fun = fun)
  return(dist)
}
dnExp <- dnExponential

dnGamma <- function(shape = 1, rate = 1) {
  min  <- 0
  max  <- qgamma(0.999, shape = shape, rate = rate)
  hpd  <- qgamma(c(0.025, 0.975), shape = shape, rate = rate)
  fun  <- function(x) dgamma(x, shape = shape, rate = rate)
  dist <- list(min = min,
               max = max,
               hpd = hpd,
               fun = fun)
  return(dist)
}










