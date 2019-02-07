#' @title read_count_data
#' @description FUNCTION_DESCRIPTION
#' @param countdata PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @param genelen PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#' @importFrom glue glue
#' @rdname read_count
#' @export
read_count <- function(countdata,
                       type,
                       genelen) {
  if (is.character(countdata)) {
    datExpr <- fread(countdata, data.table = FALSE) %>% as.matrix()
  } else {
    datExpr <- countdata %>% as.matrix()
  }

  message(glue("number of genes in raw count matrix {nrow(datExpr)}"))
  message(glue("number of cells in raw count matrix {ncol(datExpr)}"))

  if (!is.integer(datExpr)) { # if the matrix has integer counts, it is probably raw data, otherwise, probably TPM
    if (length(genelen) != nrow(datExpr)) {
      stop("number of genes in 'genelen' and count matrix do not match!")
    }
    # multiply each gene count by its gene length
    datExpr <- sweep(
      x = datExpr,
      MARGIN = 1,
      STATS = genelen,
      FUN = "*"
    )
  }

  totalCounts_by_cell <- colSums(datExpr)

  totalCounts_by_cell[totalCounts_by_cell == 0] <- 1
  datExpr <- sweep(
    x = datExpr,
    MARGIN = 2,
    STATS = 10^6 / totalCounts_by_cell,
    FUN = "*"
  )
  if (min(datExpr) < 0) {
    stop("Smallest read count cannot be negative!")
  }
  count_lnorm <- log10(datExpr + 1.01)
  return(count_lnorm)
}
