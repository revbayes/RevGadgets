#' get MAP
#'
#' Calculates the Maximum a Posteriori estimate for the trace of a
#' quantitative variable
#'
#' Uses the SANN method of the optim() function to approximate the MAP estimate
#'
#' @param var (numeric vector; no default) Vector of the samples from the
#' trace of a quantitative variable
#'
#' @return the MAP estimate
#'
#' @seealso \link[stats]{optim}
#'
#' @examples
#'
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_GTR.p", package="RevGadgets")
#' trace <- readTrace(paths = file)
#' MAP <- getMAP(trace[[1]]$"pi[1]")
#'
#' @export

getMAP <- function(var) {
  d <- density(var)
  f <- approxfun(d$x, d$y)
  op <- stats::optim(
    par = mean(var),
    fn = f,
    method = "SANN",
    control = list(fnscale = -1)
  )
  return(op$par)
}
