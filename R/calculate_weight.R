#' @title calculate_weight
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param paramt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom dplyr bind_cols
#' @importFrom stats dgamma dnorm
#' @rdname calculate_weight
#' @export
calculate_weight <- function(x, paramt) {
    pz1 <- paramt[1] * dgamma(x, 
                              shape = paramt[2], 
                              rate = paramt[3])
    pz2 <- (1 - paramt[1]) * dnorm(x, 
                                   mean = paramt[4], 
                                   sd = paramt[5])
    pz <- pz1 / (pz1 + pz2)
    pz[pz1 == 0] <- 0
    return(cbind(pz, 1 - pz))
  }
