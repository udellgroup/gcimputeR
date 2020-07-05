#' Return a complete Z matrix based on SVD initialization
#'
#' @description  The SVD initialization is based on a complete Z matrix with missing entries filled in by zero. Useful for getting initial estimates of copula parameters
#' @param Z_meanimp Initial complete Z with mean (i.e. zero) imputation
#' @param rank the number of latent factors
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @return A complte Z matrix
#' @export
impute_init = function(Z_meanimp, rank, r_lower, r_upper){
  obj = svd(Z_meanimp, nv=rank)
  Z = obj$u[,1:rank] %*% diag(obj$d[1:rank]) %*% t(obj$v)
  #ind = which(Z>r_upper | Z<r_lower)
  ind = union(which(Z>r_upper), which(Z<r_lower))
  Z[ind] = Z_meanimp[ind]
  Z
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
#' @export
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
#' @export
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
quad = function(A,x){
  as.numeric(t(x) %*% A %*% x)
}

#' Impute missing entries for the Z matrix
#'
#' @description  Using the returned quantities from EM algorithm, impute missing entries for the Z matrix
#' @param Z The returned \code{Zobs} matrix from EM algorithm
#' @param S Returned from EM algorithm
#' @param W The returned estimate for \code{W} matrix
#' @return The imputed complete \code{Z} matrix
#' @export
imputeZ_mixedgc_ppca = function(Z, W, S){
  n = dim(Z)[1]
  p = dim(Z)[2]
  obj = svd(W)
  U = obj$u
  d = obj$d
  rank = length(d)

  for (i in 1:n){
    index_m = which(is.na(Z[i,]))
    Z[i,index_m] =  matrix(U[index_m,], ncol = rank) %*% matrix(S[i,], ncol = 1)
  }
  Z
}
