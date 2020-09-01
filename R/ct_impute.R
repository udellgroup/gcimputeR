#' confidence intervals for Gaussian copula imputation
#'
#' @description  Construct confidence intervals for imputation from Gaussian copula, only applied to real-valued/continuous matrices.
#' @param X The original incomplete observation
#' @param fit The returned object from impute_mixedgc_ppca
#' @param alpha The significance level
#' @return A list containing two matrices \code{upper} and \code{lower} with the same size of \code{X}
#' \describe{
#'   \item{\code{upper}}{Upper bound of the constructed interval at missing location, \code{NA} as observed location.}
#'   \item{\code{lower}}{Lower bound of the constructed interval at missing location, \code{NA} as observed location.}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Matrix Completion with Quantified Uncertainty through Low Rank Gaussian Copula. arXiv preprint arXiv:2006.10829.
#' @export
ct_impute = function(X, fit, alpha = 0.95){
  n = dim(X)[1]
  p = dim(X)[2]
  # Fitted model parameters
  obj = svd(fit$W)
  U = obj$u
  d = obj$d
  rank = length(d)
  sig = fit$sigma
  # imputed values of Z
  Zimp = fit$Zimp

  upper = array(NA, dim = c(n,p))
  lower = array(NA, dim = c(n,p))
  var = array(NA, dim = c(n,p))
  margin = qnorm(1-(1-alpha)/2)
  # compute conditional variance
  for (i in 1:n){
    index_m = which(is.na(X[i,]))
    index_o = which(!is.na(X[i,]))
    Uiobs = matrix(U[index_o,], ncol = rank)
    Uimis = matrix(U[index_m,], ncol = rank)

    dUmis = solve(sig * diag(d^{-2}) + t(Uiobs) %*% Uiobs, t(Uimis))
    # compute variance i.e. diagonal elements of conditional covariance matrix
    for (l in 1:length(index_m)){
      j = index_m[l]
      du = dUmis[,l]
      var_ij = sig + sig * sum(du * U[j,])
      var[i,j] = var_ij
      upper[i,j] = Zimp[i,j] + margin * sqrt(var_ij)
      lower[i,j] = Zimp[i,j] - margin * sqrt(var_ij)
    }
  }

  for (j in 1:p){
    miss_ind = which(is.na(X[,j]))
    upper[miss_ind,j] = quantile(X[,j], pnorm(upper[miss_ind,j]), na.rm = TRUE)
    lower[miss_ind,j] = quantile(X[,j], pnorm(lower[miss_ind,j]), na.rm = TRUE)
  }

  return(list(upper = upper, lower = lower))
}
