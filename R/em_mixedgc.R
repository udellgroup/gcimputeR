#' EM algorithm to fit Gaussian copula
#'
#' @description  fit the Gaussian copula model from incomplete mixed data
#' @param Z_continuous Continuous columns
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param start Initial value of copula correlation matrix. Default is \code{NULL}.
#' @param maxit maximum number of iterations
#' @param eps Convergence threshold
#' @param stop.relative Use relative Frobeinus error when \code{TRUE} and maximum elementise error when \code{FALSE} by default
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
em_mixedgc = function(Z_continuous, r_lower, r_upper, start =NULL, maxit=100, eps=1e-3, stop.relative = TRUE){
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
    est_iter = em_mixedgc_iter(Z, r_lower, r_upper, rep(0,p), R)
    Z = est_iter$Zobs
    R1 = cov2cor(est_iter$sigma)
    loglik = c(loglik, est_iter$loglik)
    if (stop.relative) err = norm(R1-R, type = 'F')/norm(R, type = 'F') else err = max(abs(R1-R))
    if (err<eps) break
    if (l > maxit){
      warning('Max iter reached in EM')
      break
    }
    R = R1
  }

  return(list(R=R, loglik=loglik, Zobs = Z, Zimp = est_iter$Zimp))
}
