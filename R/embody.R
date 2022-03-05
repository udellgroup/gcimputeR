#' Compute multivariate truncated normal mean and var
#'
#' @description  iteratively compute the conditional mean and var at observed ordinal entries
#' @param z An incomplete row
#' @param lower Lower boundary of truncated intervals for ordinal columns
#' @param upper Upper boundary of truncated intervals for ordinal columns
#' @param obs_indices Boolean vector where \code{TRUE} indicates observed entries.
#' @param ord_obs_indices Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param ord_in_obs Boolean vector where \code{TRUE} indicates ordinal entries.
#' @param obs_in_ord Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param f_sigma_oo_inv_z A function computing The matrix-vector product Sigma_{obs, obs}^{-1} * z_{obs}
#' @param sigma_oo_inv_diag The diagonal of Sigma_{obs, obs}^{-1}
#' @param n_update The number of full cycle iterations to conduct
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{The updated mean at observed entries }
#'   \item{\code{var}}{The log-likelihood achieved during iteration.}
#' }
#' @export
update_z_row_ord <- function(z, lower, upper,
                             obs_indices,
                             ord_obs_indices, ord_in_obs, obs_in_ord,
                             f_sigma_oo_inv_z,
                             sigma_oo_inv_diag,
                             n_update=1){
  p = length(z)
  num_ord = length(lower)
  var_ordinal = numeric(p)

  # when there is an observed ordinal to be imputed and another observed dimension, impute this ordinal
  # The update here only requires the observed variables, but needs to use the relative location of
  # ordinal variables in all observed variables.
  if (sum(obs_indices)>1 & any(ord_obs_indices)){
    ord_obs_iter = which(ord_obs_indices)
    ord_in_obs_iter = which(ord_in_obs)
    obs_in_ord_iter = which(obs_in_ord)

    for (i in 1:n_update){
      sigma_oo_inv_z = f_sigma_oo_inv_z(z[obs_indices])
      for (index in seq_along(ord_obs_iter)){
        # j is the location in the p-dim coordinate
        # j_in_obs is the location of j in the obs-dim coordinate
        # j_in_ord is the location of j in the ord-dim coordinate
        j_in_ord = obs_in_ord_iter[index]
        j_in_obs = ord_in_obs_iter[index]
        j = ord_obs_iter[index]

        cond_var_j = 1/sigma_oo_inv_diag[j_in_obs]
        cond_std_j = sqrt(cond_var_j)
        cond_mean_j = z[j] - cond_var_j*sigma_oo_inv_z[j_in_obs]

        new_mean_j = mean_tnorm(cond_mean_j, cond_std_j, lower[j_in_ord], upper[j_in_ord])
        new_var_j = var_tnorm(cond_mean_j, cond_std_j, lower[j_in_ord], upper[j_in_ord])

        if (is.finite(new_mean_j)) z[j] = new_mean_j
        if (is.finite(new_var_j)) var_ordinal[j] = new_var_j
      }
    }

  }

  list('mean'=z, 'var'=var_ordinal)
}

#' Compute multivariate truncated normal mean and var
#'
#' @description  iteratively compute the conditional mean and var at observed ordinal entries
#' @param z An incomplete row
#' @param lower Lower boundary of truncated intervals for ordinal columns
#' @param upper Upper boundary of truncated intervals for ordinal columns
#' @param obs_indices Boolean vector where \code{TRUE} indicates observed entries.
#' @param ord_obs_indices Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param ord_in_obs Boolean vector where \code{TRUE} indicates ordinal entries.
#' @param obs_in_ord Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param sigma_oo Sigma_{obs, obs}
#' @param method Explicit derivation `Explicit` or sampling method `TruncatedNormal`
#' @param n_sample Number of MC samples to use
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{The updated mean at observed entries }
#'   \item{\code{var}}{The log-likelihood achieved during iteration.}
#' }
#' @export
est_z_row_ord <- function(z, lower, upper,
                          obs_indices, ord_obs_indices, obs_in_ord, ord_in_obs, sigma_oo,
                          method='Explicit', n_sample=5000){
  # need to determine the conditional mean and cov of observed ordinal given observed continuous
  z_obs = z[obs_indices]
  p = length(z)
  if (sum(obs_indices)>1 & any(ord_obs_indices)){
    cont_in_obs = !ord_in_obs
    sigma_ord_ord = sigma_oo[ord_in_obs, ord_in_obs, drop=FALSE]

    if (any(cont_in_obs)){
      sigma_cont_ord = sigma_oo[cont_in_obs, ord_in_obs, drop=FALSE]
      sigma_cont_cont = sigma_oo[cont_in_obs, cont_in_obs, drop=FALSE]

      tot_m = cbind(z_obs[cont_in_obs], sigma_cont_ord)
      sol_m = solve(sigma_cont_cont, tot_m)

      cond_mean = t(sigma_cont_ord) %*% sol_m[,1]
      cond_cov = sigma_ord_ord - t(sigma_cont_ord) %*% sol_m[,-1,drop=FALSE]
    }else{
      cond_mean = numeric(sum(ord_obs_indices))
      cond_cov = sigma_ord_ord
    }

    if (!isSymmetric(cond_cov)){
      diff = max(abs(cond_cov - t(cond_cov)))
      if (diff > 1e-4) print(paste('max diff:', diff))
      cond_cov = (cond_cov + t(cond_cov))/2
    }

    cond_mean = c(cond_mean)
    lower_use = lower[obs_in_ord]
    upper_use = upper[obs_in_ord]
    lmu = length(cond_mean)
    ll = length(lower_use)
    lu = length(upper_use)
    if (lmu != ncol(cond_cov) | lmu != ll | lmu != lu){
      stop('inconsistent shapes')
    }
    out_ = get_trunc_2dmoments(c(cond_mean), cond_cov,
                               lower_use, upper_use,
                               method=method, n_sample=n_sample)
    z_new = z
    z_new[ord_obs_indices] = out_$mean
    cov_all = matrix(0, p, p)
    cov_all[ord_obs_indices, ord_obs_indices] = out_$cov
    out = list('mean' = z_new, 'cov'=cov_all, 'var'=NULL)
  }else{
    out = list('mean' = z, 'var' = numeric(p))
  }

  # return: list with names 'mean' and 'cov' or NULL
  out
}


#' First two moments of truncated multivariate normal
#'
#' @description  compute multivariate truncated mean vector and covariance matrix
#' @param mean Normal mean vector
#' @param cov Normal covariance matrix
#' @param lower Lower boundary of truncated intervals
#' @param upper Upper boundary of truncated intervals
#' @param method Method to use
#' @param n_sample Number of samples to use for sampling methods
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{Mean vector}
#'   \item{\code{cov}}{Covariance matrix}
#' }
#' @export
get_trunc_2dmoments <- function(mean, cov, lower, upper, method='Explicit', n_sample=5000){
  p = length(mean)
  if (p==1){
    out = moments_truncnorm(c(mean), c(cov), c(lower), c(upper))
    return(list(mean = out$mean, cov = matrix(out$var, 1, 1)))
  }
  switch (method,
          'TruncatedNormal' = {
            z = TruncatedNormal::rtmvnorm(n = n_sample, mu = mean, sigma = cov, lb = lower, ub = upper)
            if (length(dim(z)) != 2 | ncol(z)!=p) stop('unexpected sample dimension')
            list('mean' = colMeans(z), 'cov' = cov(z))
          },
          'Explicit' = {
            r = tmvtnorm::mtmvnorm(mean = mean, sigma = cov, lower = lower, upper = upper)
            names(r) <- c('mean', 'cov')
            r
          },
          stop("Invalid `method` vlaue")
  )
}
