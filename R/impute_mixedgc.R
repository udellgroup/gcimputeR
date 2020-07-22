#' @import stats
NULL
#' @importFrom MASS mvrnorm
NULL


#' Impute mixed type data with Gaussian copula
#'
#' @description  Impute the missing values of a continuous and ordinal mixed dataset through Gaussian copula model.
#' @param X A matrix or data.frame with missing values. Observed entry of \code{X} should either be numerical value or numerical ordinal level. Make sure there is no empty row nor character level in \code{X}.
#' @param maxit Maximum number of iterations
#' @param eps Convergence threshold
#' @param stop.relative The algorithm stops when the fitted copula correlation matrix converges. The convergence is measureed using either relative Frobenius error when \code{TRUE} or maximum elementise error when \code{FALSE} by default
#' @param nlevels A column which has larger number of unique values than \code{nlevels} will be classfied as continuous, otherwise ordinal.
#' @details Impute the missing entries of a continuous and ordinal mixed data by fitting a Gaussian copula model to the data. The algorithm first scales original observation \code{X} to copula observation \code{Z} whose marginals are all standard normal. For continuous columns of \code{Z}, each observed entry records a value. While for ordinal columns of \code{Z}, each observed entry records an interval. The second step is to estimate the Gaussian copula correlation matrix using observed information of \code{Z}. With estimated correlation matrix, the third step imputes missing entries of \code{Z}. The last step imputes the missing entries of \code{X} based on the imputed entries of \code{Z}.
#' @return A list containing:
#' \describe{
#'   \item{\code{Ximp}}{Imputed data matrix}
#'   \item{\code{R}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration. This value approximates the true objective function we want to maximize, which is hard to evaluate. Monotonically increasing \code{loglik} sequence indicates good fit}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2019). Missing Value Imputation for Mixed Data Through Gaussian Copula. arXiv preprint arXiv:1910.12845.
#' @export
#' @examples
#' # Simulate Data
#' library(MASS)
#' Sigma = matrix(c(1,.8,.8,.8,1,.64,.8,.64,1), nrow = 3)
#' Z = mvrnorm(n=2000, mu = rep(0,3), Sigma = Sigma)
#' X = Z # X has the same first column as Z, as standard normal marginal
#' X[,2] = continuous2ordinal(Z[,2], k = 2) # binary marginal
#' X[,3] = continuous2ordinal(Z[,3], k = 5) # ordinal marginal with 5 levels
#'
#' # Mask 10% Data in each column
#' Xobs = X
#' Xobs[1:100,1] = NA
#' Xobs[101:200,2] = NA
#' Xobs[201:300,3] = NA
#'
#' # Fit Gaussian copula
#' fit = impute_mixedgc(Xobs)
#'
#' # Compute imputation Error
#' cal_mae_scaled(xhat = fit$Ximp, xobs = Xobs, xtrue = X)
#'

impute_mixedgc = function(X, maxit=100, eps=1e-3, stop.relative = TRUE, nlevels = 20){
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.numeric(as.matrix(X))
  dim(X) = c(n,p)

  # Do not allow empty row
  if (any(apply(X, 1, function(x){sum(!is.na(x))}) == 0)) stop("remove empty row")
  # Do not allow column with only one level
  if (any(apply(X, 2, function(x){length(unique(x[!is.na(x)]))}) <= 1)) stop('remove column with only 0 or 1 unique value')

  d_index = which(apply(X, 2, function(x) {length(unique(x)) <= nlevels}))
  k = length(d_index)
  c_index = setdiff(1:p, d_index)
  if (k > 0) {
    r = range_transform(matrix(X[, d_index], nrow = n), type = "ordinal")
    r_lower = r$r_lower
    r_upper = r$r_upper
    rm(r)
  }
  else {
    r_lower = NULL
    r_upper = NULL
  }
  if(k<p){
    Z_continuous = range_transform(matrix(X[,c_index], nrow = n),type = "continuous")$r_val
  }else{
    Z_continuous = NULL
  }
  fit_em = em_mixedgc(Z_continuous = Z_continuous, r_lower = r_lower,
                      r_upper = r_upper, maxit = maxit, eps = eps, stop.relative = stop.relative)
  R = fit_em$R
  # Impute X using Imputed Z
  Xnew.p = Ximp_transform(Z = fit_em$Zimp, X = X[, c(d_index,c_index)], d_index = d_index)
  # Back to original permuation
  R1 = R
  Xnew = Xnew.p
  if ((k > 0) & (k < p)) {
    R1[d_index, d_index] = R[1:k, 1:k]
    Xnew[, d_index] = Xnew.p[, 1:k]
    R1[d_index, c_index] = R[1:k, (k + 1):p]
    R1[c_index, d_index] = R[(k + 1):p, 1:k]
    R1[c_index, c_index] = R[(k + 1):p, (k + 1):p]
    Xnew[, c_index] = Xnew.p[, (k + 1):p]
  }
  return(list(Ximp = Xnew, R = R1, loglik = fit_em$loglik))
}

#' Impute mixed type data with low rank Gaussian copula
#'
#' @description  Impute the missing values of a continuous and ordinal mixed dataset through low rank Gaussian copula model.
#' @param X A matrix or data.frame with missing values. Observed entry of \code{X} should either be numerical value or numerical ordinal level. Make sure there is no empty row nor character level in \code{X}.
#' @param rank the number of latent factors
#' @param maxit Maximum number of iterations
#' @param eps Convergence threshold
#' @param verbose If true, output each iteration's detail
#' @param early.stop If true, stop the algorithm when the likelihood change is small.
#' @param nlevels A column which has larger number of unique values than \code{nlevels} will be classfied as continuous, otherwise ordinal.
#' @details Impute the missing entries of a continuous and ordinal mixed data by fitting a low rank Gaussian copula model to the data. The algorithm first scales original observation \code{X} to copula observation \code{Z} whose marginals are all standard normal. For continuous columns of \code{Z}, each observed entry records a value. While for ordinal columns of \code{Z}, each observed entry records an interval. The second step is to estimate the Gaussian copula correlation matrix using observed information of \code{Z}. With estimated correlation matrix, the third step imputes missing entries of \code{Z}. The last step imputes the missing entries of \code{X} based on the imputed entries of \code{Z}.
#' @return A list containing:
#' \describe{
#'   \item{\code{Ximp}}{Imputed data matrix}
#'   \item{\code{W}}{Fitted latent low rank subspace matrix}
#'   \item{\code{sigma}}{Fitted noise variance}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration. This value approximates the true objective function we want to maximize, which is hard to evaluate. Monotonically increasing \code{loglik} sequence indicates good fit}
#'   \item{\code{Zimp}}{The imputed Z matrix. On observed ordinal entries, the entry is the corresponding estimated conditional mean. Useful for constructing confidence intervals.}
#'   \item{\code{C}}{The conditional variance corresponding to the observed Z matrix. Useful for quantifying imputation uncertainty.}
#'   \item{\code{cutoffs}}{The estimated cutoffs for ordinal dimensions. Useful for quantifying imputation uncertainty.}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Matrix Completion with Quantified Uncertainty through Low Rank Gaussian Copula. arXiv preprint arXiv:2006.10829.
#' @export

impute_mixedgc_ppca = function(X, rank, maxit=50, eps=1e-6, nlevels = 20,verbose = FALSE, early.stop = FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.numeric(as.matrix(X))
  dim(X) = c(n,p)

  # Do not allow empty row
  if (any(apply(X, 1, function(x){sum(!is.na(x))}) == 0)) stop('remove empty row')
  # Do not allow column with only one level
  if (any(apply(X, 2, function(x){length(unique(x[!is.na(x)]))}) <= 1)) stop('remove column with only 0 or 1 unique value')

  # find ordinal dimensions
  d_index = which(apply(X, 2, function(x){length(unique(x))<=nlevels}))
  k = length(d_index)
  # find continuous dimensions
  c_index = setdiff(1:p, d_index)
  # marginal transformation
  if (k>0){
    r = range_transform(matrix(X[,d_index], nrow = n), type = 'ordinal') # matrix of latent scalars corresponding to measured values
    r_lower = r$r_lower
    r_upper = r$r_upper
    cutoffs = list(k)
    for (j in 1:k){
      c = unique(r$r_lower[,j])
      c = c[!is.na(c)] # remove NA
      cutoffs[[j]] = c[is.finite(c)] # remove -Inf
    }
    rm(r)
  }else{
    r_lower = NULL
    r_upper = NULL
    cutoffs = NULL
  }

  if(k<p){
    Z_continuous = range_transform(matrix(X[,c_index], nrow = n), type = 'continuous')$r_val
  }else{
    Z_continuous = NULL
  }

  # EM: estimate correlation matrix
  start = NULL
  #start = list(sig = 0, W = scale.corr(W =  matrix(rnorm(p*rank), ncol = rank), sig = 0)$W)
  fit_em = em_mixedgc_ppca(rank = rank,
                             Z_continuous =  Z_continuous, r_lower = r_lower, r_upper=r_upper,
                             maxit = maxit, eps = eps, early.stop = early.stop,
                             start = start,verbose = verbose)

  W = fit_em$W
  sigma = fit_em$sigma

  # impute missing Z
  Zimp = imputeZ_mixedgc_ppca(Z = fit_em$Zobs, S = fit_em$S, W = W)

  # Impute X using Imputed Z
  Xnew.p = Ximp_transform(Z = Zimp, X = X[,c(d_index,c_index)], d_index = d_index)

  # Back to original permuation
  W1 = W
  Xnew = Xnew.p
  if((k>0) & (k<p)){
    W1[d_index,] = W[1:k,]
    Xnew[,d_index] = Xnew.p[,1:k]
    W1[c_index,] = W[(k+1):p,]
    Xnew[,c_index] = Xnew.p[,(k+1):p]
  }
  #return(list(Ximp=Xnew, W = W1, sigma = sigma, loglik=fit_em$loglik, Zimp = Zimp, C = fit_em$C, cutoffs = cutoffs))
  return(list(Ximp=Xnew, W = W1, sigma = sigma, loglik=fit_em$loglik, Zimp = Zimp, C = fit_em$C, cutoffs = cutoffs))
}
