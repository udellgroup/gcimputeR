#' Each iteration of EM algorithm fit
#'
#' @description  fit the Gaussian copula model with incomplete mixed type data
#' @param Z Incomplete Z matrix
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param mu Mean vector. Default setting forces it to be \code{0}.
#' @param sigma Covariance matrix.
#' @param n_update The number of updates to conduct when `method='Iterative'`
#' @param trunc_method Method for evaluating truncated normal moments
#' @param n_sample Number of samples to use for sampling methods for evaluating truncated normal moments
#' @return A list containing fitted copula correlation matrix, the likelihood(objective function), Z matrix with updated ordinal entries and a complete imputed Z matrix.
#' \describe{
#'   \item{\code{R}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Zobs}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Missing Value Imputation for Mixed Data Through Gaussian Copula. KDD 2020.
#' @export
em_mixedgc_iter = function(Z, r_lower, r_upper, mu = rep(0, ncol(Z)), sigma, n_update=1,
                           trunc_method='Iterative', n_sample=5000){
  # input: "Z" matrix with missing values
  n = dim(Z)[1]
  p = dim(Z)[2]
  if (is.null(r_lower)) k = 0 else k = dim(r_lower)[2]
  negloglik = 0
  Zimp = matrix(NA, n, p)
  C = matrix(0, nrow = p, ncol = p)

  for (i in 1:n){
    missing_indices = is.na(Z[i,])
    obs_indices = !missing_indices

    # TODO: ord_indices should be input later after canceling the ordering permutation
    ord_indices = (1:p) <= k
    ord_obs_indices = ord_indices & obs_indices
    ord_in_obs = ord_obs_indices[obs_indices]
    obs_in_ord = ord_obs_indices[ord_indices]

    sigma_oo = sigma[obs_indices,obs_indices, drop=FALSE]
    sigma_om = sigma[obs_indices,missing_indices, drop=FALSE]
    sigma_mm = sigma[missing_indices,missing_indices, drop=FALSE]

    n_obs = sum(obs_indices)
    if (any(missing_indices)){
      ans = solve(sigma_oo, cbind(diag(n_obs), sigma_om))
      sigma_oo_inv = ans[,1:n_obs, drop=FALSE]
      J_mo = t(ans[,-(1:n_obs), drop=FALSE])
    }else{
      sigma_oo_inv = solve(sigma_oo)
      J_mo = NULL
    }

    # ORDINAL DIMENSION
    # Only implement when there is at least one observed ordinal dimension(to be imputed)
    # and another observed dimension (ordinal or continuous) used as information to impute
    switch (trunc_method,
            'Iterative' = {
              f_sigma_oo_inv_z <- function(zobs){sigma_oo_inv %*% zobs}
              out <- update_z_row_ord(Z[i,], r_lower[i,], r_upper[i,],
                                      obs_indices = obs_indices,
                                      ord_obs_indices = ord_obs_indices,
                                      ord_in_obs = ord_in_obs,
                                      obs_in_ord = obs_in_ord,
                                      f_sigma_oo_inv_z = f_sigma_oo_inv_z,
                                      sigma_oo_inv_diag = diag(sigma_oo_inv),
                                      n_update = n_update
              )
            },
            'Explicit' = {
              out <- est_z_row_ord(Z[i,], r_lower[i,], r_upper[i,],
                                   obs_indices = obs_indices,
                                   ord_obs_indices = ord_obs_indices,
                                   ord_in_obs = ord_in_obs,
                                   obs_in_ord = obs_in_ord,
                                   sigma_oo = sigma_oo,
                                   method = 'Explicit'
                                   )
            },
            'Sampling_TN' = {
              out <- est_z_row_ord(Z[i,], r_lower[i,], r_upper[i,],
                                   obs_indices = obs_indices,
                                   ord_obs_indices = ord_obs_indices,
                                   ord_in_obs = ord_in_obs,
                                   obs_in_ord = obs_in_ord,
                                   sigma_oo = sigma_oo,
                                   method = 'TruncatedNormal',
                                   n_sample = n_sample
              )
            },
            stop('invalid trunc_method')
    )

    Z[i,] = out$mean
    if (any(is.na(Z[i,obs_indices]))) stop('invalid Zobs')
    if (is.null(out$cov)){
      if (is.null(out$var)) stop('wrong return from trunc ordinal update')
      cov_ordinal = diag(out$var)
    }else{
      cov_ordinal = out$cov
    }
    C = C + cov_ordinal


    # Record pseudo-loglikelihood
    z_obs = Z[i,obs_indices,drop=FALSE]
    mu_obs = mu[obs_indices]
    negloglik = negloglik + c(determinant(sigma_oo)$modulus) + (z_obs-mu_obs) %*% (sigma_oo_inv %*% t(z_obs-mu_obs))
    negloglik = negloglik + p*log(2*pi)

    Zimp[i,obs_indices] = z_obs
    # MISSING DIMENSION
    if (any(missing_indices)) {
      # mean expectation: also imputation
      Zimp[i,missing_indices] = mu[missing_indices] + J_mo %*% t(z_obs - mu_obs)
      if (any(is.na(Zimp[i,]))) stop('invalid imputation')

      # variance expectation
      C[missing_indices, missing_indices] = C[missing_indices, missing_indices] + sigma_mm - J_mo %*% sigma_om
      if (sum(diag(cov_ordinal))>0){
        #n_obs_ord = sum(ord_obs_indices)
        #cov_missing_obs_ord = J_mo[,ord_in_obs, drop=FALSE] %*% diag(diag(cov_ordinal)[ord_obs_indices], n_obs_ord, n_obs_ord)
        cov_missing_obs_ord = J_mo[,ord_in_obs, drop=FALSE] %*% cov_ordinal[ord_obs_indices, ord_obs_indices, drop=FALSE]
        C[missing_indices,ord_obs_indices] = C[missing_indices,ord_obs_indices] + cov_missing_obs_ord
        C[ord_obs_indices,missing_indices] = C[ord_obs_indices,missing_indices] + t(cov_missing_obs_ord)
        C[missing_indices,missing_indices] = C[missing_indices,missing_indices] + cov_missing_obs_ord %*% t(J_mo[,ord_in_obs,drop=FALSE])
      }
    }
  }

  # variance expecation on discrete*discrete
  C = C/n
  #mu = apply(Zimp, 2, mean)
  sigma = cov(Zimp) + C
  if (any(is.na(sigma))) stop('invalid sigma')
  loglik = -negloglik/(2*n)
  return(list(sigma=sigma, loglik=loglik, Zimp=Zimp, Zobs=Z))
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
    obs_indices = which(!is.na(Z[i,]))
    d_in_o = which(obs_indices <= k)
    index_d = obs_indices[d_in_o]

    #
    zi_obs = matrix(Z[i,obs_indices], ncol = 1)
    Ui_obs = matrix(U[obs_indices,],ncol = rank)
    UU_obs = t(Ui_obs) %*% Ui_obs

    # used in both ordinal and factor block
    ans = solve(sigma * diag(d^{-2}) + UU_obs, cbind(diag(rank), t(Ui_obs)))
    AU = ans[,-(1:rank)]
    Ai = ans[,1:rank]
    A[i,,] = Ai


    # Only implement when there is at least one observed ordinal dimension(to be imputed)
    # and another observed dimension (ordinal or continuous) used as information to impute
    if (length(obs_indices)>1 & length(index_d) >=1){
      mu = (zi_obs - Ui_obs %*% (AU %*% zi_obs))/sigma # length is observed length < p

      for (j in index_d){
        j_in_o = which(obs_indices == j)
        ind_j = setdiff(obs_indices, j)

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
    zi_obs = matrix(Z[i,obs_indices], ncol = 1)
    si = matrix(AU %*% zi_obs, ncol = 1)
    S[i,] = si
    tSC[i,,] = si %*% t(si)
    #SC[i,,] =  AU %*% diag(C[i,obs_indices], nrow = length(obs_indices))  %*% t(AU)
    SC[i,,] = AU %*% (C[i,obs_indices] * t(AU))

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
    obs_indices = which(!is.na(Z[i,]))
    Wi_obs = matrix(Wnew[obs_indices,], ncol = rank)
    zi_obs = matrix(Z[i,obs_indices], ncol = 1)
    wsi = matrix(Wnew[obs_indices,], ncol = rank) %*% matrix(S[i,], ncol=1) # O_i * k
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

