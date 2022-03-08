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
#'   \item{\code{corr}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The average log-likelihood achieved during iteration.}
#'   \item{\code{Z}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#'   \item{\code{C}}{The conditional co-variance due to missingness, i.e. E(z_i z_j|x_i x_j)}
#'   \item{\code{var_ordinal}}{The conditional variance due to truncation, i.e. E(z|a < z < b)}
#' }
em_mixedgc_iter = function(Z, r_lower, r_upper, mu = rep(0, ncol(Z)), sigma, n_update=1,
                           trunc_method='Iterative', n_sample=5000){
  n = nrow(Z)
  p = ncol(Z)

  out = list('Z'=Z, 'Zimp'=Z, 'loglik'=0, 'C'=matrix(0,p,p), 'var_ordinal'=matrix(0,n,p))

  for (i in 1:n){
    row_out = latent_operation_row(Z[i,], r_lower[i,], r_upper[i,],
                               mu=mu, sigma=sigma,
                               trunc_method = trunc_method,
                               n_update=n_update, n_sample=n_sample
                               )
    for (name in c('Z', 'Zimp', 'var_ordinal')){
      out[[name]][i,] = row_out[[name]]
    }
    for (name in c('loglik', 'C')){
      out[[name]] = out[[name]] + row_out[[name]]/n
    }
  }

  out[['corr']] = cov(out[['Zimp']]) + out[['C']]
  if (any(is.na(out[['corr']]))) stop('invalid correlation')
  out
}

#' Operation for each latent row
#'
#' @description  operation for each latent row
#' @param z Incomplete Z matrix
#' @param lower Lower boundary of truncated intervals for ordinal columns
#' @param upper Upper boundary of truncated intervals for ordinal columns
#' @param mu Mean vector. Default setting forces it to be \code{0}.
#' @param sigma Covariance matrix.
#' @param n_update The number of updates to conduct when `method='Iterative'`
#' @param trunc_method Method for evaluating truncated normal moments
#' @param n_sample Number of samples to use for sampling methods for evaluating truncated normal moments
#' @return A list containing fitted copula correlation matrix, the likelihood(objective function), Z matrix with updated ordinal entries and a complete imputed Z matrix.
#' \describe{
#'   \item{\code{loglik}}{The average log-likelihood achieved during iteration.}
#'   \item{\code{Z}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#'   \item{\code{C}}{The conditional co-variance due to missingness, i.e. E(z_i z_j|x_i x_j)}
#'   \item{\code{var_ordinal}}{The conditional variance due to truncation, i.e. E(z|a < z < b)}
#' }
latent_operation_row <- function(z, lower, upper, mu, sigma, n_update=1, trunc_method='Iterative', n_sample=5000){
  out_return = list()

  missing_indices = is.na(z)
  obs_indices = !missing_indices
  # TODO: ord_indices should be input later after canceling the ordering permutation
  p = ncol(sigma)
  if (is.null(lower)) k = 0 else k = length(lower)
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
            f_sigma_oo_inv_z <- function(zz){sigma_oo_inv %*% zz}
            out <- update_z_row_ord(z, lower, upper,
                                    obs_indices = obs_indices,
                                    ord_obs_indices = ord_obs_indices,
                                    ord_in_obs = ord_in_obs,
                                    obs_in_ord = obs_in_ord,
                                    f_sigma_oo_inv_z = f_sigma_oo_inv_z,
                                    sigma_oo_inv_diag = diag(sigma_oo_inv),
                                    n_update = n_update)
          },
          'Explicit' = {
            out <- est_z_row_ord(z, lower, upper,
                                 obs_indices = obs_indices,
                                 ord_obs_indices = ord_obs_indices,
                                 ord_in_obs = ord_in_obs,
                                 obs_in_ord = obs_in_ord,
                                 sigma_oo = sigma_oo,
                                 method = 'Explicit'
            )
          },
          'Sampling_TN' = {
            out <- est_z_row_ord(z, lower, upper,
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

  # truncated moments
  z = out$mean
  if (any(is.na(z[obs_indices]))) stop('invalid Zobs')
  if (is.null(out$cov)){
    if (is.null(out$var)) stop('wrong return from trunc ordinal update')
    cov_ordinal = diag(out$var)
  }else{
    cov_ordinal = out$cov
  }
  out_return[['Z']] = z
  out_return[['var_ordinal']] = diag(cov_ordinal)

  # loglik
  z_obs = z[obs_indices]
  mu_obs = mu[obs_indices]
  negloglik = c(determinant(sigma_oo)$modulus) + sum((z_obs-mu_obs) * c(sigma_oo_inv %*% (z_obs-mu_obs)))
  negloglik = negloglik + p*log(2*pi)
  loglik = -negloglik/2
  out_return[['loglik']] = loglik

  # imputation
  zimp = z
  if (any(missing_indices)) {
    zimp[missing_indices] = mu[missing_indices] + J_mo %*% (z_obs - mu_obs)
  }
  if (any(is.na(zimp))) stop('invalid imputation')
  out_return[['Zimp']] = zimp

  C = cov_ordinal
  # MISSING DIMENSION
  if (any(missing_indices)) {
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
  out_return[['C']] = C

  out_return
}
