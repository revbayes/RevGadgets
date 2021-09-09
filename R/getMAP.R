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
#' \donttest{
#' # download the example dataset to working directory
#' url <-
#'   "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.log"
#' dest_path <- "primates_cytb_GTR.log"
#' download.file(url, dest_path)
#'
#' # to run on your own data, change this to the path to your data file
#' file <- dest_path
#'
#' trace <- readTrace(paths = file)
#' MAP <- getMAP(trace[[1]]$"pi[1]")
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path)
#' }
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
