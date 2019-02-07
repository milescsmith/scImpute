#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datExpr PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: 'count'
#' @param drop_thre PARAM_DESCRIPTION, Default: 0.5
#' @param Kcluster PARAM_DESCRIPTION, Default: NULL
#' @param labels PARAM_DESCRIPTION, Default: NULL
#' @param genelen PARAM_DESCRIPTION, Default: NULL
#' @param plan PARAM_DESCRIPTION, Default: "multiprocess"
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname scimpute
#' @export
scimpute <- function(datExpr, 
                     type = "count",
                     drop_thre = 0.5, 
                     Kcluster = NULL, 
                     labels = NULL, 
                     genelen = NULL, 
                     plan = "multiprocess") {
    
    if (is.null(labels) & is.null(Kcluster)) {
      stop("'Kcluster' must be specified when no labels are provided.")
    }
    if (!(type %in% c("count", "TPM"))) {
      stop("expression values must be either 'count' or 'TPM'.")
    }
    if (type == "TPM" & is.null(genelen)) {
      stop("For TPM, 'genelen' must be specified.")
    }

    message("reading in raw count matrix ...")
    count_lnorm <- read_count(countdata = datExpr, 
                              type = type, 
                              genelen = genelen)

    message("starting imputation...")
    if (is.null(labels)) {
      res_imp <- imputation_model8(count = count_lnorm, 
                                   labeled = FALSE,
                                   point = log10(1.01), 
                                   drop_thre = drop_thre,
                                   Kcluster = Kcluster)
    } else {
      if (length(labels) != ncol(count_lnorm)) {
        stop("number of cells does not match number of labels !")
      }
      res_imp <- imputation_wlabel_model8(count = count_lnorm,
                                          cell_labels = labels, 
                                          point = log10(1.01),
                                          drop_thre = drop_thre,
                                          Kcluster = NULL)
    }
    
    res_imp$count_imp <- 10^res_imp$count_imp - 1.01
    rownames(res_imp$count_imp) <- rownames(count_lnorm)
    colnames(res_imp$count_imp) <- colnames(count_lnorm)
    return(res_imp)
  }
