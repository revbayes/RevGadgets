#' @export
plotPrior <- function(distribution, col) {

  # check for a variable name
  if ( grepl("~", distribution) ) {
    title        <- distribution
    split        <- strsplit(distribution, "~")[[1]]
    param_name   <- split[1]
    param_name   <- gsub(" ", "", param_name)
    distribution <- tail(split, n=1)
  } else {
    title      <- paste0("x ~ ", distribution)
    param_name <- "x"
  }

  # transform Rev vectors to R vectors
  string <- distribution
  distribution <- gsub("v(", "c(", distribution, fixed=TRUE)
  distribution <- gsub("[", "c(", distribution, fixed=TRUE)
  distribution <- gsub("]", ")", distribution, fixed=TRUE)

  # get the distribution
  distn <- .getPrior(distribution)

  if ( distn$type == "continuous"  ) {

    # plot continuous distributions this way

    # get colors
    if ( missing(col) ) {
      col <- colFun(length(distn$fun))
    } else {
      # TODO: check the length of the provided color vector
    }

    # make the plots
    pp <- ggplot2::ggplot(NULL)
    for(i in 1:length(distn$fun)) {
      if (length(distn$fun) > 1) {
        pp <- pp + ggplot2::stat_function(fun = distn$fun[[i]], color = col[i], aes(color = paste0(param_name,"[",i,"]")), n=1001)
      } else {
        pp <- pp + ggplot2::stat_function(fun = distn$fun[[i]], color = col[i], aes(color = param_name), n=1001)
      }
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
    pp <- pp + ggplot2::xlim(distn$min, distn$max)

  } else if ( distn$type == "discrete" ) {

    # plot discrete distributions this way

    # get colors
    if ( missing(col) ) {
      col <- colFun(1)
    } else {
      # TODO: check the length of the provided color vector
    }

    # make the plots
    pp <- ggplot2::ggplot(distn$probs, aes(x = x, y = y, fill = is_credible)) + ggplot2::geom_bar(stat="identity", color = col)

    # fill the credible set
    pp <- pp + ggplot2::scale_fill_manual(values = alpha(c("TRUE" = col, "FALSE" = NA), 0.5))

    # scale the x axis
    pp <- pp + ggplot2::scale_x_continuous(breaks = pretty(distn$probs$x))

  } else {
    stop("invalid distribution type")
  }

  # style
  pp <- pp + ggplot2::theme_bw()

  # axis labels
  if ( length(distn$fun) > 1 ) {
    pp <- pp + ggplot2::labs(x = paste0(param_name, "[i]"), y = paste0("prior probability of ", param_name, "[i]"), title = title)
  } else {
    pp <- pp + ggplot2::labs(x = param_name, y = paste0("prior probability of ", param_name), title = title)
  }

  # title and legend, if necessary
  if ( distn$type == "continuous" ) {
    pp <- pp + ggplot2::theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  } else {
    pp <- pp + ggplot2::theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  }

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

ln <- log

###########################
# chi-square distribution #
###########################

dnChisq <- function(df) {

  min  <- 0
  max  <- qchisq(0.999, df = df)
  fun  <- function(x) dchisq(x, df = df)

  # create the HPDs
  hpd <- HDInterval::hdi(qchisq, df=df)

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

  return(dist)

}


############################
# exponential distribution #
############################

dnExponential <- function(rate = 1, offset = 0) {

  min  <- offset + 0
  max  <- offset + qexp(0.999, rate = rate)
  fun  <- function(x) dexp(x - offset, rate = rate)

  # create the HPDs
  hpd <- HDInterval::hdi(qexp, rate=rate)

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

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

  # create the HPDs
  hpd <- HDInterval::hdi(qgamma, shape = shape, rate = rate)

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

  return(dist)

}

##########################
# lognormal distribution #
##########################

dnLognormal <- function(mean = 0, sd = 1, offset = 0) {

  min  <- offset + 0
  max  <- offset + qlnorm(0.999, mean = mean, sd = sd)
  fun  <- function(x) dlnorm(x - offset, mean = mean, sd = sd)

  # create the HPDs
  hpd <- HDInterval::hdi(qlnorm, mean = mean, sd = sd) + offset

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

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

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

  return(dist)

}

dnNorm <- dnNormal

#######################
# cauchy distribution #
#######################

dnCauchy <- function(location, scale) {

  min  <- qcauchy(0.01, location = location, scale = scale)
  max  <- qcauchy(0.99, location = location, scale = scale)
  hpd  <- qcauchy(c(0.025, 0.975), location = location, scale = scale)
  fun  <- function(x) dcauchy(x, location = location, scale = scale)

  dist <- list(min  = min,
               max  = max,
               hpd  = list(hpd),
               fun  = list(fun),
               type = "continuous")

  return(dist)

}

###############################
# normal mixture distribution #
###############################

dnBimodalNormal <- function(mean1, mean2, sd1, sd2, probability) {

  # create the function
  fun <- function(x) probability * dnorm(x, mean = mean1, sd = sd1) + (1 - probability) * dnorm(x, mean = mean2, sd = sd2)

  # compute the quantiles
  cdf <- function(x) probability * pnorm(x, mean = mean1, sd = sd1) + (1 - probability) * pnorm(x, mean = mean2, sd = sd2)
  quantiles <- HDInterval::inverseCDF(c(0.001, 0.999), cdf)
  min <- quantiles[1]
  max <- quantiles[2]

  # create a density object
  x <- seq(min, max, length.out = 1001)
  y <- fun(x)
  dens <- list(x = x, y = y)
  class(dens) <- "density"

  # create the HPDs
  hdi <- HDInterval::hdi(dens, allowSplit=TRUE)

  # make the object
  dist <- list(min  = min,
               max  = max,
               hpd  = lapply(1:nrow(hdi), function(x) hdi[x,]),
               fun  = list(fun),
               type = "continuous")

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
  hdi <- HDInterval::hdi(dens, allowSplit=TRUE)

  # make the object
  dist <- list(min  = bottom,
               max  = top,
               hpd  = lapply(1:nrow(hdi), function(x) hdi[x,]),
               fun  = list(fun),
               type = "continuous")
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

  # create the HPDs
  hdi <- HDInterval::hdi(qbeta, shape1 = alpha, shape2 = beta)

  # make the object
  dist <- list(min  = min,
               max  = max,
               hpd  = list(hdi),
               fun  = list(fun),
               type = "continuous")

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

    # create the HPDs
    HDInterval::hdi(qbeta, shape1 = this_alpha, shape2 = this_beta)

  })

  # make the object
  dist <- list(min  = min,
               max  = max,
               hpd  = hpd,
               fun  = fun,
               type = "continuous")
  return(dist)

}


#########################
# binomial distribution #
#########################

.getDiscreteCredibleSet <- function(probs, cutoff = 0.95) {
  sort_probs <- sort(probs, decreasing=TRUE)
  cutoff <- min(which(cumsum(sort_probs) > cutoff))
  if (cutoff > length(sort_probs)) {
    cutoff <- length(sort_probs)
  }
  is_credible <- probs >= sort_probs[cutoff]
  return(is_credible)
}

dnBinomial <- function(p, n) {

  # the limits
  min  <- qbinom(0.001, size = n, prob = p)
  max  <- qbinom(0.999, size = n, prob = p)

  # the range
  range <- min:max
  probs <- dbinom(range, size = n, prob = p)

  # in credible set
  is_credible <- .getDiscreteCredibleSet(probs)

  # make the object
  dist <- list(min   = min,
               max   = max,
               probs = data.frame(x = range, y = probs, is_credible = is_credible),
               type  = "discrete")

  return(dist)

}


########################
# poisson distribution #
########################

dnPoisson <- function(lambda) {

  # the limits
  min  <- qpois(0.001, lambda = lambda)
  max  <- qpois(0.999, lambda = lambda)

  # the range
  range <- min:max
  probs <- dpois(range, lambda = lambda)

  # in credible set
  is_credible <- .getDiscreteCredibleSet(probs)

  # make the object
  dist <- list(min   = min,
               max   = max,
               probs = data.frame(x = range, y = probs, is_credible = is_credible),
               type  = "discrete")

  return(dist)

}

dnPois <- dnPoisson


##########################
# geometric distribution #
##########################

dnGeometric <- function(p) {

  # the limits
  min  <- 0
  max  <- qgeom(0.999, prob = p)

  # the range
  range <- min:max
  probs <- dgeom(range, prob = p)

  # in credible set
  is_credible <- .getDiscreteCredibleSet(probs)

  # make the object
  dist <- list(min   = min,
               max   = max,
               probs = data.frame(x = range, y = probs, is_credible = is_credible),
               type  = "discrete")

  return(dist)

}

dnGeom <- dnGeometric

##################################
# negative binomial distribution #
##################################

dnNbinomial <- function(r, p) {

  # the limits
  min  <- qnbinom(0.001, size = r, prob = p)
  max  <- qnbinom(0.999, size = r, prob = p)

  # the range
  range <- min:max
  probs <- dnbinom(range, size = r, prob = p)

  # in credible set
  is_credible <- .getDiscreteCredibleSet(probs)

  # make the object
  dist <- list(min   = min,
               max   = max,
               probs = data.frame(x = range, y = probs, is_credible = is_credible),
               type  = "discrete")

  return(dist)

}

dnNbinom <- dnNbinomial













