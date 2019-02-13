#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param pars PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom stats dgamma dnorm
#' @rdname dmix
#' @export
dmix <- function(x, pars) {
    pars[1] * dgamma(x, shape = pars[2], 
                     rate = pars[3]) + 
    (1 - pars[1]) * 
    dnorm(x, 
          mean = pars[4], 
          sd = pars[5])
  }
