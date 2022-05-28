#' Mask data
#'
#' @description Masks \code{mask_fraction} entries of \code{X} completely at random
#' @param X data to be masked
#' @param mask_fraction Fraction of observed entries to be masked
#' @param seed Seed for mask
#' @param silence_rate The random masking is done while ensuring there is no empty row, which may lead to smaller than specified \code{mask_fraction}. If the difference is more than \code{silence_rate}, a message will be printed
#' @param mask_cols If not \code{NULL}, the masking only happens in columns \code{mask_cols}
#' @return masked data
#' @export
mask_MCAR = function(X, mask_fraction, seed=NULL,silence_rate=0.01, mask_cols = NULL){
  if (is.null(dim(X))) return(mask_MCAR_vec(X, mask_fraction, seed))
  if (!is.null(seed)) set.seed(seed)
  X = to_numeric_matrix(X)
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
  X_masked
}

mask_MCAR_vec <- function(X, mask_fraction, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  obs_coors = which(!is.na(X))
  num  = length(obs_coors)
  mask_indices = sample(1:num, num*mask_fraction)
  mask_coors = obs_coors[mask_indices]
  X_masked = X
  X_masked[mask_coors] = NA
  X_masked
}

# val_ratio is the ratio of validation to training
split_mask_val_test <- function(X_mask, X, val_ratio = 0.5, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  X = to_numeric_matrix(X)
  if (val_ratio == 0){
    list(train = X_mask, test = X)
  }else{
    list(train = mask_MCAR(X_mask, mask_fraction = val_ratio), validation = X_mask, test = X)
  }
}

split_mask_val_test_ <- function(X_mask, X, val_ratio = 0.5, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  if (val_ratio == 0) return(list(test = X))
  X = to_numeric_matrix(X)
  mask_coors = which(!is.na(X) & is.na(X_mask))
  num  = length(mask_coors)
  index_val = sample(1:num, num*val_ratio)
  index_test = setdiff(1:num, index_val)

  r = list(validation = X_mask, test = X_mask)
  loc = list(validation = index_val, test = index_test)
  for (name in names(loc)){
    coors = mask_coors[loc[[name]]]
    r[[name]][coors] = X[coors]
  }
  r
}
