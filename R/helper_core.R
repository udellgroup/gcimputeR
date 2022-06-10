observed <- function(x) x[!is.na(x)]

sum_list_len <- function(s) sum(purrr::map_int(s, length))


#' DataFrame to Matrix
#'
#' @description Safely turn a numerical data frame to a numerical matrix
#' @param X a data frame
#' @export
#' @keywords internal
to_numeric_matrix <- function(X){
  Xnew = as.numeric(as.matrix(X))
  dim(Xnew) = dim(X)
  Xnew
}

index_int_to_logi <- function(index, l){
  out = logical(l)
  out[index] = TRUE
  out
}

regularize_corr <- function(sigma, corr_min_eigen = 1e-3, verbose = FALSE, prefix = ''){
  o_eigen = eigen(sigma)
  if (min(o_eigen$values)<corr_min_eigen){
    values = o_eigen$values
    values[values<corr_min_eigen] = corr_min_eigen
    sigma = cov2cor(o_eigen$vectors %*% diag(values) %*% t(o_eigen$vectors))
    if (verbose) print(paste0(prefix, 'small eigenvalue in the copula correlation'))
  }
  sigma
}


factor_to_num <- function(x) as.numeric(levels(x))[x]


to_nearest_ord <- function(x, ords=NULL, xobs=NULL){
  if (is.null(ords)) ords = unique(xobs[!is.na(xobs)])
  ords[which.min(abs(x - ords))]
}

to_nearest_ord_vec <- function(x, xobs){
  ords = unique(xobs[!is.na(xobs)])
  n = length(x)
  xnew = numeric(n)
  for (i in seq_along(x)) xnew[i] = to_nearest_ord(x[i], ords)
  xnew
}


#' Return a complete Z matrix based on SVD initialization
#'
#' @description  The SVD initialization is based on a complete Z matrix with missing entries filled in by zero. Useful for getting initial estimates of copula parameters
#' @param Z_meanimp Initial complete Z with mean (i.e. zero) imputation
#' @param rank the number of latent factors
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @return A complte Z matrix
#' @export
#' @keywords internal
impute_init = function(Z, rank, Lower, Upper, d_index){
  Z_imp = Z
  Z_imp[is.na(Z)] = 0
  obj = svd(Z_imp, nv=rank)
  Z_imp = obj$u[,1:rank] %*% diag(obj$d[1:rank]) %*% t(obj$v)

  p = ncol(Z)
  j_in_ord = 1
  for (j in 1:p){
    index_o = which(!is.na(Z[,j]))
    if (d_index[j]){
      index = Z_imp[index_o,j] < Lower[index_o,j_in_ord] | Z_imp[index_o,j] > Upper[index_o,j_in_ord]
      Z_imp[index_o[index],j] = Z[index_o[index],j]
      j_in_ord = j_in_ord + 1
    }else{
      Z_imp[index_o,j] = Z[index_o,j]
    }
  }

  Z_imp
}

#' Adjust estiamted correlation parameters
#'
#' @description  adjust estiamted correlation parameters, W and sigma, to satisfy unit diagonal constraints.
#' @param W The latent low rank subspace matrix (initial estimate)
#' @param sigma The noise variance (initial estimate)
#' @return A list containing
#' \describe{
#'   \item{\code{W}}{Adjusted latent low rank subspace matrix}
#'   \item{\code{sigma}}{Adjusted noise variance}
#' }
#' @export
#' @keywords internal
scale_corr = function(W, sigma){
  p = dim(W)[1]
  tr = apply(W, 1, function(x){sum(x^2)})
  sigma = mean(1/(tr+sigma)) * sigma
  for (j in 1:p) W[j,] = sqrt(1-sigma) * W[j,]/sqrt(tr[j])
  return(list(W = W, sigma = sigma))
}

#' Sum a 3d matrix along one axis
#'
#' @description  For a 3d matrix M with dimension n*k*k, sum over dimension (2,3) with scale c, at entries selected by index in dimension 1
#' @param M A 3d matrix with dimension n*k*k
#' @param c A vector of length n
#' @param index a subset of \{1,2,...,n\}
#' @return A k*k matrix
#' @keywords internal
sum_3d_scale = function(M, c, index){
  apply(M, c(2,3), function(x){sum(x[index] * c[index])})
}

#' Sum a 2d matrix along one axis
#'
#' @description  For a 2d matrix M with dimension n*k, sum over dimension 2 with scale c, at entries selected by index in dimension 1
#' @param M A 2d matrix with dimension n*k
#' @param c A vector of length n
#' @param index a subset of \{1,2,...,n\}
#' @return A k*1 matrix
#' @keywords internal
sum_2d_scale = function(M, c, index){
  v = apply(M, 2, function(x){sum(x[index] * c[index])})
  matrix(v,ncol = 1)
}

#' Quadratic form
#'
#' @description  Compute the quadratic form x^tAx
#' @param x A p*1 vector
#' @param A A p*p matrix
#' @return A scalar
#' @export
#' @keywords internal
quad = function(A,x){
  as.numeric(t(x) %*% A %*% x)
}

#' Impute missing entries for the Z matrix
#'
#' @description  Using the returned quantities from EM algorithm, impute missing entries for the Z matrix
#' @param Z The returned \code{Zobs} matrix from EM algorithm
#' @param W The returned estimate for \code{W} matrix
#' @return The imputed complete \code{Z} matrix
#' @export
#' @keywords internal
imputeZ_mixedgc_ppca = function(Z, W, sigma){
  n = dim(Z)[1]
  p = dim(Z)[2]
  obj = svd(W)
  U = obj$u
  d = obj$d
  rank = length(d)

  for (i in 1:n){
    index_m = is.na(Z[i,])
    index_o = !index_m
    Ui_obs = U[index_o,,drop=FALSE]
    zi_obs = matrix(Z[i, index_o], ncol = 1)
    UU_obs = t(Ui_obs) %*% Ui_obs
    si = solve(sigma * diag(d^{-2}) + UU_obs, t(Ui_obs) %*% zi_obs)
    Z[i,index_m] =  matrix(U[index_m,], ncol = rank) %*% matrix(si, ncol = 1)
  }
  Z
}


#' Conditional normal mean and cov
#' @description Compute the conditional normal mean and cov
#' @param mu normal mean
#' @param cov normal covariance
#' @param z data observation
#' @param index_o observed dimensions
#' @param index_m missing dimensions
#' @param drop Whether to drop the dimension of computed conditional mean
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{conditional mean}
#'   \item{\code{cov}}{conditional covariance }
#' }
#' @keywords internal
get_cond_dist <- function(z, mu, cov, index_o, index_m=NULL,drop=TRUE){
  test_logical = TRUE
  if (test_logical){
    if (!is.logical(index_o)) stop('use logical index')
    if (!any(index_o)) print('empty set to conditional on')
    if (is.null(index_m)) index_m = !index_o
    if (sum(index_o)+sum(index_m)>ncol(cov)) stop('inconsistent cov, index_o and index_m')
  }
  if (is.null(index_m)) stop('input index_m')
  if (is.null(dim(z))) z = matrix(z, 1)
  if (ncol(z) != length(mu)) stop('inconsistent z and mu')
  # dim |o| by n_sample, using the centered version
  if (any(index_o)){
    zo_c = t(z[,index_o,drop=FALSE]) - mu[index_o]
    if (any(is.na(zo_c))) stop('missing found in centered zo')
    cov_oo = cov[index_o, index_o, drop=FALSE]
    cov_om = cov[index_o, index_m, drop=FALSE]
    # |o| * |m|
    ans = solve(cov_oo, cov_om)
    # sample * |m|
    cond_mean = t(mu[index_m] + t(ans) %*% zo_c)
    cond_cov = cov[index_m,index_m, drop=FALSE] - t(cov_om) %*% ans
  }else{
    mu_m = mu[index_m]
    cond_mean = matrix(mu_m, nrow(z), length(mu_m), byrow=TRUE)
    cond_cov = cov[index_m,index_m, drop=FALSE]
  }

  if (drop) cond_mean = drop(cond_mean)
  list('mean'=cond_mean, 'cov'=cond_cov)
}
