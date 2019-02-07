#' @title fn
#' @description root-finding equation
#' @param alpha PARAM_DESCRIPTION
#' @param target PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname fn
#' @export
fn <- function(alpha, target) {
  log(alpha) - digamma(alpha) - target
}


#' @title update_gmm_pars
#' @description update parameters in gamma distribution
#' @param x PARAM_DESCRIPTION
#' @param wt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname update_gmm_pars
#' @export
#' @importFrom stats uniroot 
update_gmm_pars <- function(x, wt) {
  tp_s <- sum(wt)
  tp_t <- sum(wt * x)
  tp_u <- sum(wt * log(x))
  tp_v <- -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0) {
    alpha <- 20
  } else {
    alpha0 <- (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20) {
      alpha <- 20
    } else {
      alpha <- uniroot(fn, c(0.9, 1.1) * alpha0,
        target = tp_v,
        extendInt = "yes"
      )$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta <- tp_s / tp_t * alpha
  return(c(alpha, beta))
}


#' @title get_mix
#' @description estimate parameters in the mixture distribution
#' @param xdata PARAM_DESCRIPTION
#' @param point PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom stats sd
#' @rdname get_mix
#' @export
get_mix <- function(xdata, point) {
  inits <- rep(0, 5)
  inits[1] <- sum(xdata == point) / length(xdata)
  if (inits[1] == 0) {
    inits[1] <- 0.01
  }
  inits[2:3] <- c(0.5, 1)
  xdata_rm <- xdata[xdata > point]
  inits[4:5] <- c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {
    inits[5] <- 0
  }
  paramt <- inits
  eps <- 10
  iter <- 0
  loglik_old <- 0

  while (eps > 0.5) {
    wt <- calculate_weight(xdata, paramt)
    paramt[1] <- sum(wt[, 1]) / nrow(wt)
    paramt[4] <- sum(wt[, 2] * xdata) / sum(wt[, 2])
    paramt[5] <- sqrt(sum(wt[, 2] * (xdata - paramt[4])^2) / sum(wt[, 2]))
    paramt[2:3] <- update_gmm_pars(x = xdata, wt = wt[, 1])

    loglik <- sum(log10(dmix(xdata, paramt)))
    eps <- (loglik - loglik_old)^2
    loglik_old <- loglik
    iter <- iter + 1
    if (iter > 100) {
      break
    }
  }
  return(paramt)
}

#' @title get_mix_parameters
#' @description FUNCTION_DESCRIPTION
#' @param count PARAM_DESCRIPTION
#' @param point PARAM_DESCRIPTION, Default: log10(1.01)
#' @param plan Default: "multiprocess"
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom furrr future_map
#' @importFrom purrr reduce
#' @importFrom dplyr bind_rows
#'
#' @rdname get_mix_parameters
#' @export
get_mix_parameters <- function(count,
                               point = log10(1.01)) {
  count <- as.matrix(count)
  null_genes <- which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
  parslist <- future_map(1:nrow(count), function(ii) {
    if (ii %% 2000 == 0) {
      print(ii)
    }
    if (ii %in% null_genes) {
      rep(NA, 5)
    } else {
      xdata <- count[ii, ]
      paramt <- try(get_mix(xdata, point), silent = TRUE)
      if (class(paramt) == "try-error") {
        paramt <- rep(NA, 5)
      }
      paramt
    }
  })
  parslist <- reduce(parslist, rbind)
  colnames(parslist) <- c("rate", "alpha", "beta", "mu", "sigma")
  return(parslist)
}
