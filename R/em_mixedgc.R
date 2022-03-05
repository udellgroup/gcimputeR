  #' EM algorithm to fit Gaussian copula
#'
#' @description  fit the Gaussian copula model from incomplete mixed data
#' @param Z_continuous Continuous columns
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param start Initial value of copula correlation matrix. Default is \code{NULL}.
#' @param maxit maximum number of iterations
#' @param eps Convergence threshold
#' @param runiter When set as a positive integer, the algorithm will run \code{runiter} iterations exactly.
#' @param trunc_method Method for evaluating truncated normal moments
#' @param n_sample Number of samples to use for sampling methods for evaluating truncated normal moments
#' @return A list containing fitted copula correlation matrix, the likelihood(objective function), Z matrix with updated ordinal entries and a complete imputed Z matrix.
#' \describe{
#'   \item{\code{R}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Zobs}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Zobs} and missing entries imputed}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2019). Missing Value Imputation for Mixed Data Through Gaussian Copula. arXiv preprint arXiv:1910.12845.
#' @export
em_mixedgc = function(Z_continuous, r_lower, r_upper,
                      start =NULL, maxit=100, eps=1e-3, runiter=0, trunc_method='Iterative', n_sample=5000){
  if (is.null(Z_continuous)){
    p = dim(r_upper)[2]
    k = p
  }
  else{
    if (is.null(r_upper)){
      p = dim(Z_continuous)[2]
      k=0
    }
    else{
      p = dim(r_upper)[2] +  dim(Z_continuous)[2]
      k = dim(r_upper)[2]
    }
  }

  # Initialize observed ordinal values
  if (k > 0){
    Z_ordinal = initZ_ordinal(r_upper = r_upper, r_lower = r_lower)
  } else Z_ordinal = NULL

  # Initialize with mean imputation for Z_continuous
  if(is.null(start)){
    Z_meanimp = cbind(Z_ordinal, Z_continuous)
    Z_meanimp[is.na(Z_meanimp)] = 0
    R = cor(Z_meanimp)
    rm(Z_meanimp)
  }else{
    R = start$R
  }

  Z = cbind(Z_ordinal, Z_continuous)
  l=0
  loglik = NULL
  repeat{
    l = l+1
    est_iter = em_mixedgc_iter(Z, r_lower, r_upper, rep(0,p), R, trunc_method = trunc_method, n_sample=n_sample)
    Z = est_iter$Zobs
    R1 = cov2cor(est_iter$sigma)
    loglik = c(loglik, est_iter$loglik)
    if (runiter==0){
      err = norm(R1-R, type = 'F')/norm(R, type = 'F')
      if (err<eps) break
      if (l > maxit){
        warning('Max iter reached in EM')
        break
      }
    }else{
      if (l>runiter) break
    }

    R = R1
  }

  return(list(R=R, loglik=loglik, Zobs = Z, Zimp = est_iter$Zimp))
}

#' EM algorithm to fit low rank Gaussian copula
#'
#' @description  fit the Gaussian copula model from incomplete mixed data
#' @param rank the number of latent factors
#' @param Z_continuous Continuous columns
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param start Initial value of copula correlation matrix. Default is \code{NULL}.
#' @param maxit maximum number of iterations
#' @param eps Convergence threshold
#' @param verbose If true, output each iteration's detail
#' @param early.stop If true, stop the algorithm when the likelihood change is small.
#' @return A list containing fitted copula parameters, the likelihood (objective function), Z matrix with updated ordinal entries and the conditional variance corresponding to the observed Z matrix.
#' \describe{
#'   \item{\code{W}}{Fitted latent low rank subspace matrix}
#'   \item{\code{sigma}}{Fitted noise variance}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Zobs}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{C}}{The conditional variance corresponding to the observed Z matrix. Useful for quantifying imputation uncertainty.}
#'   \item{\code{S}}{Required quantity to impute the Z matrix}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Matrix Completion with Quantified Uncertainty through Low Rank Gaussian Copula. arXiv preprint arXiv:2006.10829.
#' @export
em_mixedgc_ppca = function(rank, Z_continuous, r_lower, r_upper, start =NULL, maxit=100, eps=1e-6,early.stop = FALSE, verbose = FALSE){
  if (is.null(Z_continuous)){
    p = dim(r_upper)[2]
    k = p
  }
  else{
    if (is.null(r_upper)){
      p = dim(Z_continuous)[2]
      k=0
    }
    else{
      p = dim(r_upper)[2] +  dim(Z_continuous)[2]
      k = dim(r_upper)[2]
    }
  }

  # Initialize observed ordinal values
  if (k > 0){
    Z_ordinal = initZ_ordinal(r_upper = r_upper, r_lower = r_lower)
  } else Z_ordinal = NULL

  # Initialize with SVD imputation for Z
  if(is.null(start)){
    Z_meanimp = cbind(Z_ordinal, Z_continuous)
    Z_meanimp[is.na(Z_meanimp)] = 0

    # approximate Z_meanimp by a low rank SVD
    # and correspondingly update Z_ordinal
    Z_meanimp = impute_init(Z_meanimp, rank, r_upper, r_lower)
    ind = which(!is.na(Z_ordinal))
    Z_ordinal[ind] = Z_meanimp[ind]
    Z = cbind(Z_ordinal, Z_continuous)
    Z_meanimp[!is.na(Z)] = Z[!is.na(Z)]

    R = cor(Z_meanimp)
    est = svd(R)
    rm(Z_meanimp, R)

    sigma = mean(est$d[-(1:rank)])
    W = est$u[,1:rank] %*% diag(sqrt(est$d[1:rank] - sigma))

    # scaling
    est = scale_corr(W,sigma)
    W = est$W
    sigma = est$sigma
    rm(est)
  }else{
    W = start$W
    sigma = start$sigma
  }

  Z = cbind(Z_ordinal, Z_continuous)
  l=0
  loglik = NULL
  #diffW = list()
  #estvar = NULL
  repeat{
    l = l+1
    est_iter = em_mixedgc_ppca_iter(Z, r_lower, r_upper, W, sigma)
    Z = est_iter$Zobs
    W1 = est_iter$W
    loglik = c(loglik, est_iter$loglik)
    #
    err = sum((W1-W)^2)/sum(W^2)
    #diffW[[l]] = c(err, grassman_dist(W1,W))
    #estvar = c(estvar, est_iter$sigma)
    if (verbose){
      print(paste('update error:', err, ' estimated var:',est_iter$sigma))
      print(paste('likelihood: ', est_iter$loglik))
    }
    if (err<eps) break
    if (l > maxit){
      warning('Max iter reached in EM')
      break
    }
    if (early.stop & l>1){
      r = abs((rev(loglik)[1] - rev(loglik)[2]))/abs(rev(loglik)[2])
      if (r < 0.01){
        if (verbose) print('early stop because changed likelihood below 1%')
        break
      }
    }
    W = W1
    sigma = est_iter$sigma
  }

  return(list(W=W, sigma = sigma, loglik=loglik, Zobs = Z, S = est_iter$S, C = est_iter$C))
}
