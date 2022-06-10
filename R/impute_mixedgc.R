#' @import stats
NULL
#' @importFrom MASS mvrnorm
NULL


optim_params <- function(){
  c(
    "@param maxit Maximum number of iterations",
    "@param eps Convergence threshold",
    "@param verbose Whether to print progress information",
    "@param runiter When set as a positive integer, the algorithm will run the specified number of iterations exactly."
  )
}


#' Gaussian copula for incomplete mixed data
#'
#' @description  Fit a Gaussian copula model from (continuous and ordinal) mixed data and impute the missing entries using the fitted model
#' @param X A matrix or data.frame with missing values. Observed entry of \code{X} should either be numerical value or numerical ordinal level. Make sure there is no empty row nor character level in \code{X}.
#' @param nlevels A column which has larger number of unique values than \code{nlevels} will be classfied as continuous, otherwise ordinal.
#' @eval optim_params()
#' @param trunc_method Method for evaluating truncated normal moments: \code{'Iterative'} or \code{'Sampling'}.
#' @param n_update The number of updates, only used when \code{trunc_method} is \code{'Iterative'}
#' @param n_sample Number of MC samples, only used when \code{trunc_method} is \code{'Sampling'}
#' @param corr If not \code{NULL}, impute missing values using \code{corr} as the copula correlation
#' @param n_MI The number of random samples to draw from the missing distribution.
#' @param ... Additional arguments for development use
#' @details Impute the missing entries of continuous and ordinal mixed data by fitting a Gaussian copula model to the data.
#' @return A list containing:
#' \describe{
#'   \item{\code{Ximp}}{Imputed data matrix}
#'   \item{\code{corr}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration. This value approximates the true objective function we want to maximize, which is hard to evaluate. Monotonically increasing \code{loglik} sequence indicates good fit}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Missing Value Imputation for Mixed Data via Gaussian Copula. KDD 2020
#' @export
#' @examples
#' # Simulate Data
#' library(MASS)
#' # Generate 15-dim mixed data and mask 10% observation
#' var_types = list('cont'=1:5, 'ord'=6:10, 'bin'=11:15)
#' X = generate_mixed_from_gc(var_types = var_types, n = 500)
#' Xmask = mask_MCAR(X, mask_fraction = 0.2)
#' # Fit Gaussian copula
#' fit = impute_mixedgc(Xmask, verbose = TRUE)
#'
#' # Compute imputation Error
#' cal_mae_scaled(xhat = fit$Ximp, xobs = Xmask, xtrue = X)
#'
impute_mixedgc = function(X, nlevels = 20,
                          trunc_method = 'Iterative', n_sample=5000, n_update=1,
                          maxit=50, eps=0.01, verbose=FALSE, runiter = 0, n_MI=0,
                          corr = NULL,...){
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.numeric(as.matrix(X))
  dim(X) = c(n,p)

  # Do not allow empty row
  if (any(apply(X, 1, function(x){sum(!is.na(x))}) == 0)) stop("remove empty row")
  # Do not allow column with only one level
  if (any(apply(X, 2, function(x){length(unique(x[!is.na(x)]))}) <= 1)) stop('remove column with only 0 or 1 unique value')
  d_index = apply(X, 2, function(x) {length(unique(x)) <= nlevels})
  c_index = !d_index

  # observed to latent
  obj_latent = observed_to_latent(X, d_index)
  Z = obj_latent$Z
  Lower = obj_latent$Lower
  Upper = obj_latent$Upper

  if (is.null(corr)){
    fit_em = em_mixedgc(Z, Lower, Upper, d_index,
                        maxit = maxit, eps = eps, runiter=runiter, verbose=verbose,
                        trunc_method = trunc_method, n_sample=n_sample, n_update=n_update,
                        ...)
    Zimp = fit_em$Zimp
    corr = fit_em$corr
    loglik = fit_em$loglik
  }else{
    out = latent_operation('fillup',
                           Z, Lower, Upper, d_index,
                           corr = corr,
                           n_update = n_update, n_sample = n_sample, trunc_method = trunc_method,
                           ...)
    Zimp = out$Zimp
    loglik = NULL
  }


  # Impute X using Imputed Z
  Ximp = Ximp_transform(Z = Zimp, X = X, d_index = d_index)
  return(list(Ximp = Ximp, corr = corr, loglik = loglik, d_index=which(d_index)))
}


#' Warapper
#'@export
#'@keywords internal
impute_LRGC <- function(...) impute_mixedgc_ppca(...)

#' Low rank Gaussian copula for incomplete mixed data
#'
#' @description  Fit a low rank Gaussian copula model from (continuous and ordinal) mixed data and impute the missing entries using the fitted model
#' @param rank The rank, i.e. number of latent factors
#' @inheritParams impute_mixedgc
#' @details Impute the missing entries of continuous and ordinal mixed data by fitting a low rank Gaussian copula (LRGC) model to the data. LRGC is a subclass of Gaussian copula: it requires the copula correlation matrix to have a low rank plus diagonal decomposition: \eqn{\Sigma = WW^\top + \sigma^2 \mathrm{I}_p} where \eqn{W\in \mathbb{R}\times {p\times k}} and \eqn{k<p}.
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
#' @references Zhao, Y., & Udell, M. (2020). Matrix Completion with Quantified Uncertainty through Low Rank Gaussian Copula. NeurIPS 2020.
#' @export
impute_LRGC = function(X, rank, nlevels = 20,
                       trunc_method = 'Iterative', n_sample=5000, n_update=1,
                       maxit=50, eps=0.01, verbose = FALSE, runiter = 0, ...){
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.numeric(as.matrix(X))
  dim(X) = c(n,p)

  # Do not allow empty row
  if (any(apply(X, 1, function(x){sum(!is.na(x))}) == 0)) stop('remove empty row')
  # Do not allow column with only one level
  if (any(apply(X, 2, function(x){length(unique(x[!is.na(x)]))}) <= 1)) stop('remove column with only 0 or 1 unique value')
  d_index = apply(X, 2, function(x) {length(unique(x)) <= nlevels})
  c_index = !d_index

  # observed to latent
  obj_latent = observed_to_latent(X, d_index)
  Z = obj_latent$Z
  Lower = obj_latent$Lower
  Upper = obj_latent$Upper

  # TODO: store cutoff
  if (any(d_index)){
    k = sum(d_index)
    cutoffs = list()
    for (j in which(d_index)){
      cuts = unique(Lower[,j])
      cuts = cuts[!is.na(cuts)]
      cutoffs[[as.character(j)]] = cuts[is.finite(cuts)]
    }
  }else cutoffs = NULL


  # EM: estimate correlation matrix
  #start = list(sig = 0, W = scale.corr(W =  matrix(rnorm(p*rank), ncol = rank), sig = 0)$W)
  fit_em = em_mixedgc_ppca(rank, Z,  Lower, Upper,
                           d_index = d_index,
                           maxit = maxit, eps = eps, runiter = runiter, verbose = verbose,
                           trunc_method = trunc_method, n_sample=n_sample, n_update=n_update,
                           ...)
  W = fit_em$W
  sigma = fit_em$sigma
  loglik = fit_em$loglik

  # impute missing Z
  Zimp = imputeZ_mixedgc_ppca(Z = fit_em$Z,  W = W, sigma = sigma)

  # Impute X using Imputed Z
  Ximp = Ximp_transform(Z = Zimp, X = X, d_index = d_index)

  list(Ximp=Ximp, W = W, sigma = sigma, loglik=loglik,
       Z=fit_em$Z, Zimp = Zimp, C = fit_em$C, cutoffs = cutoffs, S=fit_em$S)
}

