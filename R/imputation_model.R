#' @title find_hv_genes
#' @description FUNCTION_DESCRIPTION
#' @param count PARAM_DESCRIPTION
#' @param I PARAM_DESCRIPTION
#' @param J PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom stats quantile
#' @importFrom furrr future_map
#' @rdname find_hv_genes
#' @export
#' @importFrom furrr future_map
find_hv_genes <- function(count, 
                          I) {
  count_nzero <- future_map(1:I, 
                            function(i) {
                              setdiff(count[i, ], log10(1.01))
                            })
  
  mu <- future_map_dbl(count_nzero, mean)
  mu[is.na(mu)] <- 0
  sd <- sapply(count_nzero, sd)
  sd[is.na(sd)] <- 0
  cv <- sd / mu
  cv[is.na(cv)] <- 0
  # sum(mu >= 1 & cv >= quantile(cv, 0.25), na.rm = TRUE)
  high_var_genes <- which(mu >= 1 & cv >= quantile(cv, 0.25))
  if (length(high_var_genes) < 500) {
    high_var_genes <- 1:I
  }
  count_hv <- count[high_var_genes, ]
  return(count_hv)
}

#' @title find_neighbors
#' @description FUNCTION_DESCRIPTION
#' @param count_hv PARAM_DESCRIPTION
#' @param labeled PARAM_DESCRIPTION
#' @param J PARAM_DESCRIPTION
#' @param Kcluster PARAM_DESCRIPTION, Default: NULL
#' @param cell_labels PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom rsvd rpca
#' @importFrom furrr future_map future_map_dbl
#' @importFrom stats quantile complete.cases
#' @importFrom kernlab specc
#' @rdname find_neighbors
#' @export
find_neighbors <- function(count_hv, 
                           labeled, 
                           J, 
                           Kcluster = NULL,
                           cell_labels = NULL) {
  if (labeled == TRUE) {
    if (class(cell_labels) == "character") {
      labels_uniq <- unique(cell_labels)
      labels_mth <- 1:length(labels_uniq)
      names(labels_mth) <- labels_uniq
      clust <- labels_mth[cell_labels]
    } else {
      clust <- cell_labels
    }
    nclust <- length(unique(clust))
    message("calculating cell distances...")
    dist_list <- future_map(1:nclust, function(ll) {
      cell_inds <- which(clust == ll)
      count_hv_sub <- count_hv[, cell_inds, drop = FALSE]
      if (length(cell_inds) < 1000) {
        npc <- dimReduc(threshold = 0.4,
                        count_mat = count_hv_sub,
                        is_labeled = labeled,
                        clusters = Kclusters, 
                        scale = FALSE)
      } else {
        npc <- dimReduc(threshold = 0.6,
                        count_mat = count_hv_sub,
                        is_labeled = labeled,
                        clusters = Kclusters, 
                        k = 1000,
                        scale = FALSE)
      }
      
      if (npc < 3) {
        npc <- 3
      }
      
      mat_pcs <- t(pca$x[, 1:npc])
      
      dist_cells_list <- future_map(1:length(cell_inds), function(id1) {
        d <- future_map_dbl(1:id1, function(id2) {
          sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
          sqrt(sse)
        })
        return(c(d, rep(0, length(cell_inds) - id1)))
      })
      dist_cells <- matrix(0, nrow = length(cell_inds), ncol = length(cell_inds))
      for (cellid in 1:length(cell_inds)) {
        dist_cells[cellid, ] <- dist_cells_list[[cellid]]
      }
      dist_cells <- dist_cells + t(dist_cells)
      return(dist_cells)
    })
    
    return(list(dist_list = dist_list, clust = clust))
  } else if (labeled == FALSE) {
    ## dimensional reduction
    message("dimensional reduction...")
    if (J < 5000) {
      mat_pcs <- dimReduc(threshold = 0.4,
                      count_mat = count_hv,
                      is_labeled = labeled,
                      clusters = Kcluster,
                      scale = FALSE)
    } else {
      mat_pcs <- dimReduc(threshold = 0.6,
                      count_mat = count_hv,
                      is_labeled = labeled,
                      clusters = Kcluster,
                      k = 1000,
                      scale = FALSE)
    }
    
    ## detect outliers
    message("calculating cell distances...")
    dist_cells_list <- future_map(1:J, function(id1) {
      d <- future_map_dbl(1:id1, function(id2) {
        sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
      })
      return(c(d, rep(0, J - id1)))})
    
    dist_cells <- matrix(0, nrow = J, ncol = J)
    
    for (cellid in 1:J) {
      dist_cells[cellid, ] <- dist_cells_list[[cellid]]}
    
    dist_cells <- dist_cells + t(dist_cells)
    
    min_dist <- future_map_dbl(1:J, function(i) {
      min(dist_cells[i, -i])})
    
    iqr <- quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
    
    outliers <- which(min_dist > 1.5 * iqr + quantile(min_dist, 0.75))
    
    ## clustering
    non_out <- setdiff(1:J, outliers)
    spec_res <- specc(t(mat_pcs[, non_out]), centers = Kcluster, kernel = "rbfdot")
    message("cluster sizes:")
    message(spec_res@size)
    nbs <- rep(NA, J)
    nbs[non_out] <- spec_res
    
    return(list(dist_cells = dist_cells, clust = nbs))
  }
}

#' @title find_va_genes
#' @description FUNCTION_DESCRIPTION
#' @param parslist PARAM_DESCRIPTION
#' @param subcount PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom stats dgamma dnorm
#' @rdname find_va_genes
#' @export
find_va_genes <- function(parslist, subcount) {
  point <- log10(1.01)
  valid_genes <- which((rowSums(subcount) > point * ncol(subcount)) &
                         complete.cases(parslist))
  if (length(valid_genes) == 0) {
    return(valid_genes)
  }
  
  # find out genes that violate assumption
  mu <- parslist[, "mu"]
  
  sgene1 <- which(mu <= log10(1 + 1.01))
  
  dcheck1 <- dgamma(mu + 1, 
                    shape = parslist[, "alpha"], 
                    rate = parslist[, "beta"])
  
  dcheck2 <- dnorm(mu + 1, 
                   mean = parslist[, "mu"], 
                   sd = parslist[, "sigma"])
  
  sgene3 <- which(dcheck1 >= dcheck2 & mu <= 1)
  sgene <- union(sgene1, sgene3)
  valid_genes <- setdiff(valid_genes, sgene)
  return(valid_genes)
}

#' @title impute_nnls
#' @description FUNCTION_DESCRIPTION
#' @param Ic PARAM_DESCRIPTION
#' @param cellid PARAM_DESCRIPTION
#' @param subcount PARAM_DESCRIPTION
#' @param droprate PARAM_DESCRIPTION
#' @param geneid_drop PARAM_DESCRIPTION
#' @param geneid_obs PARAM_DESCRIPTION
#' @param nbs PARAM_DESCRIPTION
#' @param distc PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname impute_nnls
#' @export
#' @importFrom penalized predict penalized
impute_nnls <- function(Ic, 
                        cellid,
                        subcount,
                        droprate,
                        geneid_drop,
                        geneid_obs,
                        nbs, 
                        distc) {
  yobs <- subcount[, cellid]
  if (length(geneid_drop) == 0 | length(geneid_drop) == Ic) {
    return(yobs)
  }
  yimpute <- rep(0, Ic)
  
  xx <- subcount[geneid_obs, nbs]
  yy <- subcount[geneid_obs, cellid]
  ximpute <- subcount[geneid_drop, nbs]
  num_thre <- 500
  if (ncol(xx) >= min(num_thre, nrow(xx))) {
    if (num_thre >= nrow(xx)) {
      new_thre <- round((2 * nrow(xx) / 3))
    } else {
      new_thre <- num_thre
    }
    filterid <- order(distc[cellid, -cellid])[1:new_thre]
    xx <- xx[, filterid, drop = FALSE]
    ximpute <- ximpute[, filterid, drop = FALSE]
  }
  
  set.seed(cellid)
  
  nnls <- penalized(yy,
                    penalized = xx, 
                    unpenalized = ~0,
                    positive = TRUE, 
                    lambda1 = 0, 
                    lambda2 = 0,
                    maxiter = 3000, 
                    trace = FALSE)
  
  ynew <- predict(nnls, 
                  penalized = ximpute, 
                  unpenalized = ~0)[, 1]
  
  yimpute[geneid_drop] <- ynew
  
  yimpute[geneid_obs] <- yobs[geneid_obs]
  
  maxobs <- future_map(subcount, 1, max)
  
  yimpute[yimpute > maxobs] <- maxobs[yimpute > maxobs]
  
  return(yimpute)
}


#' @title imputation_model8
#' @description FUNCTION_DESCRIPTION
#' @param count PARAM_DESCRIPTION
#' @param labeled PARAM_DESCRIPTION
#' @param point PARAM_DESCRIPTION
#' @param drop_threshold PARAM_DESCRIPTION, Default: 0.5
#' @param Kcluster PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname imputation_model8
#' @export
imputation_model8 <- function(count,
                              labeled = FALSE,
                              point,
                              drop_threshold = 0.5,
                              Kcluster = 10) {
  return(imputation_model(count = count, 
                          labeled = labeled, 
                          point = point, 
                          drop_threshold = drop_threshold, 
                          Kcluster = Kcluster))
}


#' @title imputation_wlabel_model8
#' @description FUNCTION_DESCRIPTION
#' @param count PARAM_DESCRIPTION
#' @param labeled PARAM_DESCRIPTION
#' @param cell_labels PARAM_DESCRIPTION, Default: NULL
#' @param point PARAM_DESCRIPTION
#' @param drop_threshold PARAM_DESCRIPTION
#' @param Kcluster PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @importFrom glue glue
#' @rdname imputation_wlabel_model8
#' @export
imputation_wlabel_model8 <- function(count, 
                                     labeled, 
                                     cell_labels = NULL, 
                                     point, 
                                     drop_threshold,
                                     Kcluster = NULL) {
  return(imputation_model(count = count, 
                          labeled = labeled,
                          point = point, 
                          drop_threshold = drop_threshold, 
                          Kcluster = Kcluster, 
                          cell_labels = cell_labels))
}


#' imputation_model
#'
#' @param count 
#' @param point 
#' @param drop_threshold 
#' @param Kcluster 
#' @param labeled 
#' @param cell_labels 
#'
#' @return
#' @export
#' @importFrom furrr future_map future_map_dbl future_map_dfc
#' @importFrom glue glue
#' @importFrom rsvd rpca
#'
#' @examples
imputation_model <- function(count,
                             point,
                             drop_threshold = 0.5,
                             Kcluster = 10,
                             labeled = FALSE,
                             cell_labels = NULL){
  if (isTRUE(labeled)){
    Kcluster = NULL
    if (!(class(cell_labels) %in% c("character", "numeric", "integer"))) {
      stop("cell_labels should be a character or integer vector!")
    }
  }
  
  count <- as.matrix(count)
  I <- nrow(count)
  J <- ncol(count)
  count_imp <- count
  
  # find highly variable genes
  count_hv <- find_hv_genes(count, I)
  
  message("searching candidate neighbors... ")
  if (Kcluster == 1) {
    clust <- rep(1, J)
    if (J < 5000) {
      npc <- dimReduc(threshold = 0.4,
                      count_mat = count_hv,
                      clusters = Kcluster,
                      is_labeled = labeled,
                      k = 1000,
                      scale = FALSE)
    } else {
      npc <- dimReduc(threshold = 0.6,
                      count_mat = count_hv,
                      clusters = Kcluster,
                      is_labeled = labeled,
                      k = 1000,
                      scale = FALSE)
    }
    
    if (npc < 3) {
      npc <- 3
    }
    mat_pcs <- t(pca$x[, 1:npc]) # columns are cells
    
    dist_cells_list <- future_map(1:J, function(id1) {
      d <- future_map_dbl(1:id1, function(id2) {
        sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
        })
      return(c(d, rep(0, J - id1)))
      })
    dist_cells <- matrix(0, nrow = J, ncol = J)
    for (cellid in 1:J) {
      dist_cells[cellid, ] <- dist_cells_list[[cellid]]
      }
    dist_cells <- dist_cells + t(dist_cells)
  } else if (isTRUE(labeled)){
    neighbors_res <- find_neighbors(count_hv = count_hv, 
                                    labeled = TRUE, 
                                    J = J,
                                    cell_labels = cell_labels)
    dist_list <- neighbors_res$dist_list
    clust <- neighbors_res$clust
  } else {  
    message("inferring cell similarities...")
    set.seed(Kcluster)
    neighbors_res <- find_neighbors(count_hv = count_hv,
                                    labeled = FALSE,
                                    J = J,
                                    Kcluster = Kcluster)
    dist_cells <- neighbors_res$dist_cells
    clust <- neighbors_res$clust
  }
  
  # mixture model
  nclust <- sum(!is.na(unique(clust)))
  
  for(cc in 1:nclust){
    print(glue("estimating dropout probability for cluster {cc}..."))
    cells <- which(clust == cc)
    parslist <- get_mix_parameters(count = count[,cells, drop = FALSE],
                                   point = log10(1.01))
    
    message("searching for valid genes...")
    valid_genes <- find_va_genes(parslist, subcount = count[, cells])
    
    if (length(cells) > 1 & length(valid_genes) > 10) {
      subcount <- count[valid_genes, cells, drop = FALSE]
      Ic <- length(valid_genes)
      Jc <- ncol(subcount)
      parslist <- parslist[valid_genes, , drop = FALSE]
      
      droprate <- t(future_map_dfc(1:Ic, function(i) {
        wt <- calculate_weight(subcount[i, ], parslist[i, ])
        return(wt[, 1])
      }))
      
      mucheck <- sweep(subcount, 
                       MARGIN = 1, 
                       parslist[, "mu"], 
                       FUN = ">")
      
      droprate[mucheck & droprate > drop_threshold] <- 0
      
      # dropouts
      setA <- future_map(1:Jc, function(cellid) {
        which(droprate[, cellid] > drop_threshold)
      })
      
      # non-dropouts
      setB <- future_map(1:Jc, function(cellid) {
        which(droprate[, cellid] <= drop_threshold)
      })
      
      # imputation
      message(glue("imputing dropout values for cluster {cc}..."))
      if (isTRUE(labeled)){
        dlist = dist_list[[cc]]
      } else {
        dlist = dist_cells[cells, cells]
      }
      subres <- future_map(.x = 1:Jc, 
                               .progress = TRUE, 
                               .f = impute_values,
                               num_cells = Jc,
                               num_genes = Ic,
                               subcount = subcount,
                               droprate = droprate,
                               dropouts = setA,
                               non_dropouts = setB,
                               dist_list = dlist)
      count_imp[valid_genes, cells] <- subres
    }
  }
  if (isTRUE(labeled)){
    outlier <- integer(0)
  } else {
    outlier <- which(is.na(clust))
  }
  
  count_imp[count_imp < point] <- point
  return(list(count_imp = count_imp, outlier = outlier))
}


#' dimReduc
#'
#' @param threshold 
#' @param count_mat 
#' @param is_labeled 
#' @param clusters 
#' @param ... Additional parameters to pass to rpca
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom rsvd rpca
dimReduc <- function(threshold, 
                     count_mat, 
                     is_labeled, 
                     clusters, 
                     ...){
  pca <- rpca(t(count_mat), ...)
  eigs <- pca$eigvals
  var_cum <- cumsum(eigs) / sum(eigs)
  if (max(var_cum) <= threshold) {
    npc <- length(var_cum)
  } else {
    npc <- which.max(var_cum > threshold)
    if (is_labeled == FALSE) {
      npc <- max(npc, clusters)
    }
  }
  
  if (npc < 3) {
    npc <- 3
  }
  
  mat_pcs <- t(pca$x[, 1:npc]) # columns are cells
  return(mat_pcs)
}

#' impute_values
#'
#' @param cellid 
#' @param num_cells 
#' @param num_genes 
#' @param subcount 
#' @param droprate 
#' @param dropouts 
#' @param dist_list 
#' @param non_dropouts 
#'
#' @return
#' @export
#'
#' @examples
impute_values <- function(cellid,
                          num_cells,
                          num_genes,
                          subcount,
                          droprate,
                          dropouts,
                          dist_list,
                          non_dropouts){
  if (cellid %% 100 == 0) {
    print(cellid)
  }
  nbs <- which(!1:num_cells %in% cellid)
  if (length(nbs) != 0) {
    y <- try(impute_nnls(num_genes,
                         cellid = cellid, 
                         subcount, 
                         droprate, 
                         dropouts[[cellid]],
                         non_dropouts[[cellid]], 
                         nbs, 
                         distc = dist_list),
             silent = TRUE)
    if (class(y) == "try-error") {
      y <- subcount[, cellid, drop = FALSE]
    }
    return(y) 
  } else {
    return(NULL) 
  }
}