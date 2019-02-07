#' @title rmix
#' @description FUNCTION_DESCRIPTION
#' @param pars PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname rmix
#' @export
#' @importFrom stats rgamma rnorm
rmix <- function(pars, n) {
    n1 <- round(n * pars[1])
    n2 <- n - n1
    x1 <- rgamma(n1, shape = pars[2], rate = pars[3])
    x2 <- rnorm(n2, mean = pars[4], sd = pars[5])
    return(c(x1, x2))
  }
