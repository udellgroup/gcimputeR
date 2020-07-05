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


#' Each iteration of EM algorithm fit (LRGC)
#'
#' @description  fit the Gaussian copula model with incomplete mixed type data
#' @param Z Incomplete Z matrix
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param W The latent low rank subspace matrix
#' @param sigma The noise variance
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
em_mixedgc_ppca_iter = function(Z, r_lower, r_upper, W, sigma){
  # input: "Z" matrix with missing values
  n = dim(Z)[1]
  p = dim(Z)[2]
  rank = dim(W)[2]
  if (is.null(r_lower)) k = 0 else k = dim(r_lower)[2]
  negloglik = 0

  # SVD on W
  w_svd = svd(W)
  U = w_svd$u
  V = w_svd$v
  d = w_svd$d

  A = array(0, dim = c(n,rank,rank))
  SC = A
  tSC = A
  S = matrix(0, nrow = n, ncol = rank)
  C = matrix(0, nrow = n, ncol = p)

  # E-step iterate over n
  for (i in 1:n){
    # indexing
    index_o = which(!is.na(Z[i,]))
    d_in_o = which(index_o <= k)
    index_d = index_o[d_in_o]

    #
    zi_obs = matrix(Z[i,index_o], ncol = 1)
    Ui_obs = matrix(U[index_o,],ncol = rank)
    UU_obs = t(Ui_obs) %*% Ui_obs

    # used in both ordinal and factor block
    ans = solve(sigma * diag(d^{-2}) + UU_obs, cbind(diag(rank), t(Ui_obs)))
    AU = ans[,-(1:rank)]
    Ai = ans[,1:rank]
    A[i,,] = Ai


    # Only implement when there is at least one observed ordinal dimension(to be imputed)
    # and another observed dimension (ordinal or continuous) used as information to impute
    if (length(index_o)>1 & length(index_d) >=1){
      mu = (zi_obs - Ui_obs %*% (AU %*% zi_obs))/sigma # length is observed length < p

      for (j in index_d){
        j_in_o = which(index_o == j)
        ind_j = setdiff(index_o, j)

        sigma_ij = sigma/(1-quad(Ai, U[j,]))
        mu_ij = Z[i,j] - mu[j_in_o] * sigma_ij

        mu_ij_new = mean_tnorm(mu_ij, sqrt(sigma_ij), a = r_lower[i,j], b = r_upper[i,j])
        sigma_ij_new = var_tnorm(mu_ij, sqrt(sigma_ij), a = r_lower[i,j], b = r_upper[i,j])

        if (is.finite(mu_ij_new)){
          Z[i,j] = mu_ij_new # UPDATE Z: observed ordinal entry
        }
        if (is.finite(sigma_ij_new)){
          C[i,j] = sigma_ij_new          #UPDATE C
        }
      }
    }

    # update observed ordinal entries
    zi_obs = matrix(Z[i,index_o], ncol = 1)
    si = matrix(AU %*% zi_obs, ncol = 1)
    S[i,] = si
    tSC[i,,] = si %*% t(si)
    #SC[i,,] =  AU %*% diag(C[i,index_o], nrow = length(index_o))  %*% t(AU)
    SC[i,,] = AU %*% (C[i,index_o] * t(AU))

    # pseudo-likelihood
    #negloglik = negloglik + log(sigma) * p + log(det(diag(rank) + diag(d) %*% UU_obs %*% diag(d)/sigma))
    negloglik = negloglik + log(sigma) * p + log(det(diag(rank) + outer(d/sigma, d) * UU_obs))
    negloglik = negloglik + sum(zi_obs^2) - c(t(zi_obs) %*% (Ui_obs %*% si))
  }


  # M-step in W iterate over p
  Wnew = W
  s = sum(C)

  AC = array(0, dim = c(p,rank,rank))
  ASC = AC
  for (j in 1:p){
    index_j = which(!is.na(Z[,j]))
    # numerator
    AC[j,,] = sum_3d_scale(A, c = C[,j], index = index_j)
    rj = sum_2d_scale(M = S, c = Z[,j], index = index_j) + AC[j,,] %*% matrix(U[j,], ncol=1)
    # denominator
    ASC[j,,] = sum_3d_scale(SC, c = rep(1,n), index = index_j) + sum_3d_scale(A, c = rep(sigma,n), index = index_j)
    Fj = sum_3d_scale(tSC, c = rep(1,n), index = index_j) +  ASC[j,,]
    # update
    Wnew[j,] = solve(Fj, rj)
  }
  rm(tSC, SC, A)
  # M-step in sigma^2
  for (i in 1:n){
    index_o = which(!is.na(Z[i,]))
    Wi_obs = matrix(Wnew[index_o,], ncol = rank)
    zi_obs = matrix(Z[i,index_o], ncol = 1)
    wsi = matrix(Wnew[index_o,], ncol = rank) %*% matrix(S[i,], ncol=1) # O_i * k
    s = s + sum((zi_obs - wsi)^2)
  }
  for (j in 1:p){
    wj = matrix(Wnew[j,], ncol = 1)
    uj = matrix(U[j,], ncol=1)
    s = s - 2*c(t(wj) %*% AC[j,,] %*% uj) + c(t(wj) %*% ASC[j,,] %*% wj)
  }
  sigma_new = s/sum(!is.na(Z))
  #Wnew = Wnew %*% diag(d) %*% t(V)
  Wnew = Wnew %*% (d * t(V))

  # scaling
  est = scale_corr(Wnew, sigma_new)

  return(list(W = est$W, sigma = est$sigma, loglik=-negloglik/2, Zobs=Z, S = S, C=C))
}

