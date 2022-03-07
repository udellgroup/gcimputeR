#' Mask data
#'
#' @description Masks \code{mask_fraction} entries of \code{X} completely at random
#' @param X data to be masked
#' @param mask_fraction Fraction of observed entries to be masked
#' @param max_try Number of attempts to conduct mask, if an empty row appears after a masking attempt
#' @param seed Seed for mask
#' @param allow_empty_row If False, allow some masked rows to be empty
#' @return masked data
#' @export
mask_MCAR = function(X, mask_fraction, max_try=50,  seed=NULL, allow_empty_row=FALSE){
  count = 0
  X_masked = X
  n = dim(X)[1]
  p = dim(X)[2]
  obs_coors = which(!is.na(X))
  num  = length(obs_coors)
  if (!is.null(seed)) set.seed(seed)
  i=0
  while(TRUE){
    mask_indices = sample(1:num, num*mask_fraction)
    mask_coors = obs_coors[mask_indices]
    X_masked[mask_coors] = NA

    empty_row = any(apply(X_masked, 1, function(x){sum(!is.na(x))}) == 0)
    if (!empty_row) break
    i=i+1
    if (i > max_try) stop('cannot produce masking without empty row')
  }

  X_masked
}
