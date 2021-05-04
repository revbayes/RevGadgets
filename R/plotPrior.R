#' @export
plotPrior <- function(distribution, col) {

  # transform Rev vectors to R vectors
  string <- distribution
  distribution <- gsub("v(", "c(", distribution, fixed=TRUE)
  distribution <- gsub("[", "c(", distribution, fixed=TRUE)
  distribution <- gsub("]", ")", distribution, fixed=TRUE)

  # get the distribution
  distn <- .getPrior(distribution)

  # get colors
  if ( missing(col) ) {
    col <- colFun(length(distn$fun))
  } else {
    # TODO: check the length of the provided color vector
  }

  # plot the distribution
  pp <- ggplot2::ggplot(data.frame())
  for(i in 1:length(distn$fun)) {
    pp <- pp + ggplot2::stat_function(fun = distn$fun[[i]], col = col[i], n=1001)
  }

  # add the HPDs
  for(i in 1:length(distn$hpd)) {
    this_interval <- distn$hpd[[i]]
    if ( length(distn$fun) > 1 ) {
      pp <- pp + ggplot2::stat_function(fun = distn$fun[[i]], xlim = this_interval, geom="area", fill=col[i], alpha=0.5, n=1001)
    } else {
      pp <- pp + ggplot2::stat_function(fun = distn$fun[[1]], xlim = this_interval, geom="area", fill=col, alpha=0.5, n=1001)
    }
  }

  # limits
  pp <- pp + ggplot2::xlim(distn$min, distn$max) + ggplot2::theme_bw()

  # style
  pp <- pp + ggplot2::theme_bw()

  # axis
  if ( length(distn$fun) > 1 ) {
    pp <- pp + ggplot2::labs(x = "x[i]", y = "prior probability of x[i]", title = string)
  } else {
    pp <- pp + ggplot2::labs(x = "x", y = "prior probability of x", title = string)
  }
  pp <- pp + ggplot2::theme(plot.title = element_text(hjust = 0.5))

  # return the plot
  return(pp)

}

.getPrior <- function(distribution) {

  # parse the distribution string
  dist_string <- parse(text = distribution)[[1]]

  # evaluate it
  dist <- eval(dist_string)

  # return the distribution object
  return(dist)

}

############################
# exponential distribution #
############################

dnExponential <- function(rate = 1, offset = 0) {

  min  <- offset + 0
  max  <- offset + qexp(0.999, rate = rate)
  fun  <- function(x) dexp(x - offset, rate = rate)

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hpd <- hdi(dens, allowSplit=TRUE)

  dist <- list(min = min,
               max = max,
               hpd = list(hpd),
               fun = list(fun))
  return(dist)

}

dnExp <- dnExponential

######################
# gamma distribution #
######################

dnGamma <- function(shape = 1, rate = 1) {

  min  <- 0
  max  <- qgamma(0.999, shape = shape, rate = rate)
  fun  <- function(x) dgamma(x, shape = shape, rate = rate)

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hpd <- hdi(dens, allowSplit=TRUE)

  dist <- list(min = min,
               max = max,
               hpd = list(hpd),
               fun = list(fun))
  return(dist)

}

##########################
# lognormal distribution #
##########################

dnLognormal <- function(mean = 0, sd = 1, offset = 0) {

  min  <- offset + 0
  max  <- offset + qlnorm(0.999, mean = mean, sd = sd)
  fun  <- function(x) dlnorm(x - offset, mean = mean, sd = sd)

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hpd <- hdi(dens, allowSplit=TRUE)

  dist <- list(min = min,
               max = max,
               hpd = list(hpd),
               fun = list(fun))
  return(dist)

}

dnLnorm <- dnLognormal

#######################
# normal distribution #
#######################

dnNormal <- function(mean = 0, sd = 1, offset = 0) {
  min  <- offset + qnorm(0.001, mean = mean, sd = sd)
  max  <- offset + qnorm(0.999, mean = mean, sd = sd)
  hpd  <- offset + qnorm(c(0.025, 0.975), mean = mean, sd = sd)
  fun  <- function(x) dnorm(x - offset, mean = mean, sd = sd)
  dist <- list(min = min,
               max = max,
               hpd = list(hpd),
               fun = list(fun))
  return(dist)
}

dnNorm <- dnNormal

###############################
# normal mixture distribution #
###############################

dnBimodalNormal <- function(mean1, mean2, sd1, sd2, probability) {

  # create the function
  fun <- function(x) probability * dnorm(x, mean = mean1, sd = sd1) + (1 - probability) * dnorm(x, mean = mean2, sd = sd2)

  # compute the quantiles
  cdf <- function(x) probability * pnorm(x, mean = mean1, sd = sd1) + (1 - probability) * pnorm(x, mean = mean2, sd = sd2)
  quantiles <- inverseCDF(c(0.001, 0.999), cdf)
  min <- quantiles[1]
  max <- quantiles[2]

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hdi <- hdi(dens, allowSplit=TRUE)

  # make the object
  dist <- list(min = min,
               max = max,
               hpd = lapply(1:nrow(hdi), function(x) hdi[x,]),
               fun = list(fun))
  return(dist)

}

########################
# soft-bounded uniform #
########################

dnSoftBoundUniformNormal <- function(min, max, sd, p) {

  # create the function
  lm <- min
  lx <- max
  s  <- sd
  lp <- p
  fun <- function(x) {
    probs <- numeric(length(x))
    probs[x < lm] <- (1 - lp) * 0.5 * dnorm(x[x < lm] - lm, sd = s)
    probs[x > lx] <- (1 - lp) * 0.5 * dnorm(x[x > lx] - lx, sd = s)
    probs[x >= lm & x <= lx] <- lp * dunif(x[x >= lm & x <= lx], lm, lx)
    return(probs)
  }

  # compute the quantiles
  bottom <- min + qnorm(0.001, 0, sd)
  top    <- max + qnorm(0.999, 0, sd)

  # create a density object
  x <- seq(bottom, top, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hdi <- hdi(dens, allowSplit=TRUE)

  # make the object
  dist <- list(min = bottom,
               max = top,
               hpd = lapply(1:nrow(hdi), function(x) hdi[x,]),
               fun = list(fun))
  return(dist)

}

#####################
# beta distribution #
#####################

dnBeta <- function(alpha, beta) {

  # range of plot
  min <- 0
  max <- 1

  # the density of function
  fun <- function(x) dbeta(x, shape1 = alpha, shape2 = beta)

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hdi <- hdi(dens, allowSplit=TRUE)

  # make the object
  dist <- list(min = min,
               max = max,
               hpd = lapply(1:nrow(hdi), function(x) hdi[x,]),
               fun = list(fun))
  return(dist)

}

##########################
# dirichlet distribution #
##########################

dnDirichlet <- function(alpha) {

  # range of plot
  min <- 0
  max <- 1

  # the sum of all alpha values
  total_alpha <- sum(alpha)

  # the density functions
  fun <- lapply(1:length(alpha), function(i) {
    this_alpha <- alpha[i]
    this_beta  <- total_alpha - this_alpha
    function(x) dbeta(x, this_alpha, this_beta)
  })

  # the HPDs
  hpd <- lapply(1:length(alpha), function(i) {

    this_alpha <- alpha[i]
    this_beta  <- total_alpha - this_alpha

    # create a density object
    x <- seq(min, max, length.out = 1001)
    y <- fun[[i]](x)
    dens <- list(x = x, y = y)
    class(dens) <- "density"

    # create the HPDs
    hdi(dens, allowSplit=TRUE)

  })

  # make the object
  dist <- list(min = min,
               max = max,
               hpd = hpd,
               fun = fun)
  return(dist)

}




















