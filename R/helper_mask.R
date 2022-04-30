#' Mask data
#'
#' @description Masks \code{mask_fraction} entries of \code{X} completely at random
#' @param X data to be masked
#' @param mask_fraction Fraction of observed entries to be masked
#' @param seed Seed for mask
#' @param allow_empty_row If False, allow some masked rows to be empty
#' @return masked data
#' @export
mask_MCAR = function(X, mask_fraction, seed=NULL,
                     allow_empty_row=FALSE, silence_rate=0.01, mask_cols = NULL){
  if (!is.null(seed)) set.seed(seed)
  count = 0
  X_masked = X
  n = dim(X)[1]
  p = dim(X)[2]
  obs_coors = which(!is.na(X))
  if (!is.null(mask_cols)){
    coor2d = arrayInd(obs_coors, .dim = c(n,p))
    in_cols = purrr::map_lgl(coor2d[,2], ~ .x %in% mask_cols)
    obs_coors = obs_coors[in_cols]
  }
  num  = length(obs_coors)
  mask_indices = sample(1:num, num*mask_fraction)
  mask_coors = obs_coors[mask_indices]
  X_masked[mask_coors] = NA
  which_empty = which(apply(X_masked, 1, function(x){sum(!is.na(x))}) == 0)
  if (!allow_empty_row){
    for (row in which_empty){
      obs_loc = which(!is.na(X[row,]))
      if (length(obs_loc)==1) index = obs_loc else index = sample(obs_loc,1)
      X_masked[row, index] = X[row, index]
    }
    if (is.null(mask_cols)) r = sum(is.na(X_masked) & !is.na(X))/sum(!is.na(X))
    else r = sum(is.na(X_masked[,mask_cols]) & !is.na(X[,mask_cols]))/sum(!is.na(X[,mask_cols]))
    if (r<mask_fraction-silence_rate){
      print(paste0('Actual masking ratio ', round(r,4)))
    }
  }
  X_masked
}

split_mask_val_test <- function(X_mask, X, val_ratio = 0.5, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  mask_coors = which(!is.na(X) & is.na(X_mask))
  num  = length(mask_coors)
  index_val = sample(1:num, num*val_ratio)
  index_test = setdiff(1:num, index_val)
  X_val = X_mask
  X_test = X_mask
  r = list(validation = X_mask, test = X_mask)
  loc = list(validation = index_val, test = index_test)
  for (name in names(loc)){
    coors = mask_coors[loc[[name]]]
    r[[name]][coors] = X[coors]
  }
  r
}


