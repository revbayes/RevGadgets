#' Color Function
#'
#' Produce default RevGadgets colors
#'
#' Produces a vector of colors from the default RevGadgets colors
#' of length given by n, maximum of 12 colors.
#'
#' @param n (integer; no default) Number of colors to return. Maximum of 12.
#'
#' @return Character vector of color hex codes.
#'
#' @examples
#'
#' my_colors <- colFun(2)
#'
#' @export
#'

colFun <- function(n) {
  if (n %% 1 != 0) {
    stop("n must be an integer")
  }
  if (n == 1) {
    return("#005ac8")
  }
  if (n == 2) {
    return(c("#005ac8", "#fa7850"))
  }
  if (n == 3) {
    return(c("#14d2dc", "#005ac8", "#fa7850"))
  }
  if (n == 4) {
    return(c("#14d2dc", "#005ac8", "#fa7850", "#aa0a3c"))
  }
  if (n == 5) {
    return(c("#14d2dc", "#005ac8", "#fa7850", "#aa0a3c",
             "#0ab45a"))
  }
  if (n == 6) {
    return(c(
      "#14d2dc",
      "#005ac8",
      "#fa7850",
      "#aa0a3c",
      "#0ab45a",
      "#006e82"
    ))
  }
  if (n == 7) {
    return(c(
      "#14d2dc",
      "#005ac8",
      "#fa7850",
      "#aa0a3c",
      "#0ab45a",
      "#006e82",
      "#fa78fa"
    ))
  }
  if (n == 8) {
    return(
      c(
        "#14d2dc",
        "#005ac8",
        "#fa7850",
        "#aa0a3c",
        "#0ab45a",
        "#006e82",
        "#fa78fa",
        "#8214a0"
      )
    )
  }
  if (n == 9) {
    return(
      c(
        "#14d2dc",
        "#005ac8",
        "#fa7850",
        "#aa0a3c",
        "#0ab45a",
        "#006e82",
        "#fa78fa",
        "#8214a0",
        "#fae6be"
      )
    )
  }
  if (n == 10) {
    return(
      c(
        "#14d2dc",
        "#005ac8",
        "#fa7850",
        "#aa0a3c",
        "#0ab45a",
        "#006e82",
        "#fa78fa",
        "#8214a0",
        "#fae6be",
        "#00a0fa"
      )
    )
  }
  if (n == 11) {
    return(
      c(
        "#14d2dc",
        "#005ac8",
        "#fa7850",
        "#aa0a3c",
        "#0ab45a",
        "#006e82",
        "#fa78fa",
        "#8214a0",
        "#fae6be",
        "#00a0fa",
        "#f0f032"
      )
    )
  }

  if (n == 12) {
    return(
      c(
        "#14d2dc",
        "#005ac8",
        "#fa7850",
        "#aa0a3c",
        "#0ab45a",
        "#006e82",
        "#fa78fa",
        "#8214a0",
        "#fae6be",
        "#00a0fa",
        "#f0f032",
        "#a0fa82"
      )
    )
  }
  if (n >= 13) {
    stop("more than 12 colors is not supported")
  }
}
