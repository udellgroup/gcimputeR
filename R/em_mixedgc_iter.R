#' Each iteration of EM algorithm fit
#'
#' @description  fit the Gaussian copula model with incomplete mixed type data
#' @param Z Incomplete Z matrix
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param mu Mean vector. Default setting forces it to be \code{0}.
#' @param sigma Covariance matrix.
#' @return A list containing fitted copula correlation matrix, the likelihood(objective function), Z matrix with updated ordinal entries and a complete imputed Z matrix.
#' \describe{
#'   \item{\code{R}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Zobs}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2019). Missing Value Imputation for Mixed Data Through Gaussian Copula. arXiv preprint arXiv:1910.12845.
#' @export
em_mixedgc_iter = function(Z, r_lower, r_upper, mu = rep(0, dim(Z)[2]), sigma){
  # input: "Z" matrix with missing values
  n = dim(Z)[1]
  p = dim(Z)[2]
  if (is.null(r_lower)) k = 0 else k = dim(r_lower)[2]
  C = matrix(0, nrow = p, ncol = p)
  negloglik = 0
  Zimp = Z

  for (i in 1:n){
    index_o = which(!is.na(Z[i,]))
    index_m = setdiff(1:p, index_o)
    d_in_o = which(index_o <= k)
    index_d = index_o[d_in_o]

    sigma11 = matrix(sigma[index_o,index_o], nrow = length(index_o))
    sigma12 = matrix(sigma[index_o,index_m], nrow = length(index_o))
    sigma22 = matrix(sigma[index_m,index_m], nrow = length(index_m))

    ans = solve(sigma11, cbind(diag(length(index_o)), sigma12))
    sigma11_inv = ans[,1:length(index_o)]
    J21 = t(matrix(ans[,-(1:length(index_o))], nrow = length(index_o)))

    var_ordinal = numeric(p)

    # ORDINAL DIMENSION
    # Only implement when there is at least one observed ordinal dimension(to be imputed)
    # and another observed dimension (ordinal or continuous) used as information to impute
    if (length(index_o)>1 & length(index_d) >=1){
      for (j in index_d){
        j_in_o = which(index_o == j)
        ind_j = setdiff(index_o, j)

        v = sigma11_inv[,j_in_o]
        sigma_ij = 1/v[j_in_o]
        mu_ij = mu[j] + c(t(v[-j_in_o]) %*% (Z[i,ind_j]-mu[ind_j])) * (-sigma_ij)

        mu_ij_new = mean_tnorm(mu_ij, sqrt(sigma_ij), a = r_lower[i,j], b = r_upper[i,j])
        var_ij_new = var_tnorm(mu_ij, sqrt(sigma_ij), a = r_lower[i,j], b = r_upper[i,j])
        if (is.finite(mu_ij_new)) Z[i,j] = mu_ij_new # mean expectation on ordinal: also imputation
        if (is.finite(var_ij_new)){
          var_ordinal[j] = var_ij_new
          C[j,j] = C[j,j] + var_ij_new # variance expectation on ordinal
        }
      }
    }

    # MISSING DIMENSION
    if (length(index_m) > 0) {
      z_obs = matrix(Z[i,index_o], ncol=1)

      # Record pseudo-loglikelihood
      negloglik = negloglik + log(abs(det(sigma11))) + t(z_obs-mu[index_o]) %*% (sigma11_inv %*% (z_obs-mu[index_o]))
      # mean expecatation: also imputation
      Zimp[i,index_m] = mu[index_m] + J21 %*% (z_obs - mu[index_o])

      # variance expecatation
      if (length(index_d) >= 1 & length(index_o) > 1 & sum(var_ordinal)>0){
        diag_var_ordinal = diag(length(index_d))
        diag(diag_var_ordinal) = var_ordinal[index_d]
        J21varAdj = matrix(J21[,d_in_o], ncol = length(index_d)) %*% diag_var_ordinal
        #
        C[index_m,index_m] = C[index_m,index_m] + (sigma22 - J21 %*% sigma12) +
          J21varAdj %*% t(matrix(J21[,d_in_o], ncol = length(index_d)))
        #
        C[index_m,index_d] = C[index_m,index_d] + J21varAdj
        #
        C[index_d,index_m] = C[index_d,index_m] + t(J21varAdj)
      } else{
        C[index_m,index_m] = C[index_m,index_m] + (sigma22 - J21 %*% sigma12)
      }
    }
  }

  # variance expecation on discrete*discrete
  C = C/n
  #mu = apply(Zimp, 2, mean)
  sigma = cov(Zimp) + C
  return(list(sigma=sigma, loglik=-negloglik, Zimp=Zimp, Zobs=Z))
}
