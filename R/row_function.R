
#' Operation to conduct in the latent space for each row
#'
#' @description  operation for each latent row
#' @param z A row of \code{Z}
#' @param lower A row of \code{Lower}
#' @param upper A row of \code{Upper}
#' @inheritParams latent_operation
#' @return A list containing
#' \describe{
#'   \item{\code{loglik}}{Available when \code{task = 'em'}. The average log-likelihood.}
#'   \item{\code{Z}}{Available when \code{task = 'em'}. Incomplete \code{Z} with updated observed ordinal entries}
#'   \item{\code{Zimp}}{Available when \code{task = 'em'} or \code{task == 'fillup'} . Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#'   \item{\code{Zimp_sample}}{Available when \code{task = 'sample'}. Multiple imputation samples.}
#'   \item{\code{C}}{Available when \code{task = 'em'}. The conditional co-variance due to missingness}
#'   \item{\code{var_ordinal}}{Available when \code{task = 'em'} or \code{task = 'fillup'}. The conditional variance due to truncation, i.e. Var(z|a < z < b)}
#' }
#' @keywords internal
latent_operation_row <- function(task,
                                 z, lower, upper,
                                 d_index, dcat_index,
                                 corr,
                                 cat_input = NULL,
                                 trunc_method='Iterative', n_sample=5000, n_update=1, n_MI=1){
  # TODO: remove var alias
  sigma = corr
  p = ncol(sigma)
  out_return = list()
  mis_indices = is.na(z)
  obs_indices = !mis_indices

  ord_indices = d_index
  ord_obs_indices = ord_indices & obs_indices
  ord_in_obs = ord_obs_indices[obs_indices]
  obs_in_ord = ord_obs_indices[ord_indices]


  sigma_oo = sigma[obs_indices,obs_indices, drop=FALSE]
  sigma_om = sigma[obs_indices,mis_indices, drop=FALSE]
  sigma_mm = sigma[mis_indices,mis_indices, drop=FALSE]

  n_obs = sum(obs_indices)
  if (any(mis_indices)){
    ans = solve(sigma_oo, cbind(diag(n_obs), sigma_om))
    sigma_oo_inv = ans[,1:n_obs, drop=FALSE]
    J_mo = t(ans[,-(1:n_obs), drop=FALSE])
  }else{
    sigma_oo_inv = solve(sigma_oo)
    J_mo = NULL
  }

  # loglik
  if (task == 'em'){
    z_obs = z[obs_indices]
    negloglik = c(determinant(sigma_oo)$modulus) + sum(z_obs * c(sigma_oo_inv %*% z_obs))
    negloglik = negloglik + p*log(2*pi)
    loglik = -negloglik/2
    out_return[['loglik']] = loglik
  }

  # preprocessing stage when categorical data exists
  if (is.null(dcat_index)) dcat_index = logical(p)
  cat_obs_indices = dcat_index & obs_indices
  cat_in_obs = cat_obs_indices[obs_indices]
  if (any(cat_obs_indices)){
    stopifnot(!is.null(cat_input))
    x_cat = cat_input[['x_cat']]
    cat_index_list = cat_input[['cat_index_list']]
    cat_obs = !is.na(x_cat)
    # A should have dim |cat_o| * |cat_o|
    A = x_to_A(x = x_cat[cat_obs], cat_index_list = cat_index_list[cat_obs],
               d_cat = sum(cat_obs_indices), old = cat_input[['old']]) #***
    # assuming A satisfying A %*% A = I
  }else A = NULL

  # ORDINAL DIMENSION
  # Only implement when there is at least one observed ordinal dimension(to be imputed)
  # and another observed dimension (ordinal or continuous) used as information to impute
  if (trunc_method == 'TruncatedNormal' | trunc_method == 'Sampling_TN') trunc_method = 'Sampling'

  "if (trunc_method == 'Iterative'){
    out_ref <- est_z_row_ord(z, lower, upper,
                             obs_indices = obs_indices,
                             ord_obs_indices = ord_obs_indices,
                             ord_in_obs = ord_in_obs,
                             obs_in_ord = obs_in_ord,
                             sigma_oo = A_sigma_tA_at_cat(sigma_oo, A, cat_index=ord_in_obs),
                             method = 'TruncatedNormal',
                             n_sample = n_sample)
    z = out_ref$mean
    if (!is.null(A)) z[ord_obs_indices] = A %*% z[ord_obs_indices]
    #out$var = diag(out_ref$cov)
  }"

  # update ordinal latent Z

  if (task %in% c('em', 'fillup')){
    switch (trunc_method,
            'Iterative' = {
              if (!is.null(A)){
                "sigma_oo_inv = A_sigma_tA_at_cat(sigma_oo_inv, A,
                                               cat_index=ord_in_obs, A_at_left = FALSE)"
                sigma_oo_inv = solve(A_sigma_tA_at_cat(sigma_oo, A, cat_index=cat_in_obs))
                z[cat_obs_indices] = A %*% z[cat_obs_indices]
              }
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
            'Sampling' = {
              out <- est_z_row_ord(z, lower, upper,
                                   obs_indices = obs_indices,
                                   ord_obs_indices = ord_obs_indices,
                                   ord_in_obs = ord_in_obs,
                                   obs_in_ord = obs_in_ord,
                                   sigma_oo = A_sigma_tA_at_cat(sigma_oo, A, cat_index=cat_in_obs),
                                   n_sample = n_sample)
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
      #if (trunc_method == 'Sampling_TN') cov_ordinal = diag(diag(cov_ordinal))
    }

    if (!is.null(A)){
      z[cat_obs_indices] = A %*% z[cat_obs_indices]
      cov_ordinal = A_sigma_tA_at_cat(cov_ordinal, A, cat_index = cat_obs_indices)
    }

    switch (task,
            'em' = {
              out_return[['Z']] = z
              out_return[['var_ordinal']] = diag(cov_ordinal)
            },
            'fillup' = {out_return[['var_ordinal']] = diag(cov_ordinal)}
    )
  }

  z_obs = z[obs_indices]
  # imputation
  zimp = z
  if (any(mis_indices)) {
    zimp[mis_indices] = J_mo %*% z_obs
  }
  if (any(is.na(zimp))) stop('invalid imputation')
  if (task %in% c('em', 'fillup')) out_return[['Zimp']] = zimp

  # sample
  if (task == 'sample'){
    # multiple imputation step 1: truncated normal sampling
    if (trunc_method == 'Sampling'){
      zsample <- sample_z_row_ord(z, lower, upper,
                                  obs_indices = obs_indices,
                                  ord_obs_indices = ord_obs_indices,
                                  ord_in_obs = ord_in_obs,
                                  obs_in_ord = obs_in_ord,
                                  sigma_oo = A_sigma_tA_at_cat(sigma_oo, A, cat_index=cat_in_obs),
                                  n_sample = n_MI)

      if (!is.null(A)){
        # zsample[,cat_obs_indices] = t(A %*% t(zsample[,cat_obs_indices,drop=FALSE]))
        zsample[,cat_obs_indices] = zsample[,cat_obs_indices,drop=FALSE] %*% t(A)
      }
    }else{
      stopifnot(n_MI == 1)
      zsample = matrix(z, 1)
    }

    # multiple imputation step 2: normal sampling
    if (any(mis_indices)) {
      mis_cond_cov = sigma_mm - J_mo %*% sigma_om
      pmis = sum(mis_indices)
      zbar = MASS::mvrnorm(n_MI, mu = rep(0, pmis), Sigma = mis_cond_cov)
      mis_cond_mean = zsample[,obs_indices,drop=FALSE] %*% t(J_mo)
      stopifnot(dim(zbar) == dim(mis_cond_mean))
      zsample[, mis_indices] = zbar + mis_cond_mean
    }

    out_return[['Zimp_sample']] = zsample
  }

  if (task == 'em'){
    C = cov_ordinal
    # MISSING DIMENSION
    if (any(mis_indices)) {
      # variance expectation
      C[mis_indices, mis_indices] = C[mis_indices, mis_indices] + sigma_mm - J_mo %*% sigma_om
      if (sum(diag(cov_ordinal))>0){
        cov_missing_obs_ord = J_mo[,ord_in_obs, drop=FALSE] %*% cov_ordinal[ord_obs_indices, ord_obs_indices, drop=FALSE]
        C[mis_indices,ord_obs_indices] = C[mis_indices,ord_obs_indices] + cov_missing_obs_ord
        C[ord_obs_indices,mis_indices] = C[ord_obs_indices,mis_indices] + t(cov_missing_obs_ord)
        C[mis_indices,mis_indices] = C[mis_indices,mis_indices] + cov_missing_obs_ord %*% t(J_mo[,ord_in_obs,drop=FALSE])
      }
    }
    out_return[['C']] = C
  }

  out_return
}

latent_operation_LRGC_row <- function(task,
                                      z, lower, upper,
                                      d_index,
                                      U,d,sigma,
                                      cat_input = NULL, dcat_index=NULL,
                                      trunc_method='Iterative', n_sample=5000, n_update=1, n_MI=1){
  out_return = list()

  p = nrow(U)
  rank = ncol(U)

  mis_indices = is.na(z)
  obs_indices = !mis_indices
  ord_indices = d_index
  ord_obs_indices = ord_indices & obs_indices
  ord_in_obs = ord_obs_indices[obs_indices]
  obs_in_ord = ord_obs_indices[ord_indices]

  #
  Ui_obs = U[obs_indices,,drop=FALSE]
  UU_obs = t(Ui_obs) %*% Ui_obs
  #
  ans = solve(sigma * diag(d^{-2}) + UU_obs, cbind(diag(rank), t(Ui_obs)))
  Ai = ans[,1:rank,drop=FALSE]
  AU = ans[,-(1:rank),drop=FALSE]

  # preprocessing stage when categorical data exists
  if (is.null(dcat_index)) dcat_index = logical(p)
  cat_obs_indices = dcat_index & obs_indices
  cat_in_obs = cat_obs_indices[obs_indices]
  if (any(cat_obs_indices)){
    stop('unimplemented')
  }else A = NULL

  stopifnot(trunc_method %in% c('Iterative', 'Sampling'))

  # update ordinal latent Z
  if (task %in% c('em', 'fillup')){
    switch (trunc_method,
            'Iterative' = {
              f_sigma_oo_inv_z = function(zobs) (zobs - (Ui_obs %*% (AU %*% zobs)))/sigma
              sigma_oo_inv_diag = (1 - quad_mul(Ai, Ui_obs))/sigma
              out <- update_z_row_ord(z, lower, upper,
                                      obs_indices = obs_indices,
                                      ord_obs_indices = ord_obs_indices,
                                      ord_in_obs = ord_in_obs,
                                      obs_in_ord = obs_in_ord,
                                      f_sigma_oo_inv_z = f_sigma_oo_inv_z,
                                      sigma_oo_inv_diag = sigma_oo_inv_diag,
                                      n_update = n_update)
            },
            'Sampling' = {
              stop('unimplemented')
            }
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

    zi_obs = z[obs_indices]
    si = matrix(AU %*% zi_obs, ncol = 1)
    var_ordinal = diag(cov_ordinal)
    switch (task,
            'em' = {
              out_return[['var_ordinal']] = var_ordinal
              out_return[['Z']] = z
              out_return[['A']] = Ai
              out_return[['S']] = si
              out_return[['SS']] =  AU %*% (var_ordinal[obs_indices] * t(AU)) + si %*% t(si)

              # pseudo-likelihood
              nobs = sum(obs_indices)
              negloglik = p*log(2*pi) + log(sigma) * (nobs-rank)
              negloglik = negloglik + c(determinant(diag(rank) * sigma + outer(d, d) * UU_obs)$modulus)
              negloglik = negloglik + sum(zi_obs^2) - (t(zi_obs) %*% (Ui_obs %*% si))/sigma
              loglik = -negloglik/2
              out_return[['loglik']] = loglik
              out_return[['zobs_norm']] = sum(zi_obs^2)
            },
            'fillup' = {
              out_return[['var_ordinal']] = var_ordinal
              zimp = z
              zimp[mis_indices] = U[mis_indices,,drop=FALSE] %*% si
              out_return[['Zimp']] = zimp
            }
    )
  }else{
    stop('unimplemented')
  }

  out_return
}





#' Compute the observed ordinal mean and var in a row
#'
#' @description  Iteratively update an estimate of the conditional mean and var at observed ordinal entries
#' @inheritParams latent_operation_row
#' @param obs_indices Boolean vector where \code{TRUE} indicates observed entries.
#' @param ord_obs_indices Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param ord_in_obs Boolean vector where \code{TRUE} indicates ordinal entries.
#' @param obs_in_ord Boolean vector where \code{TRUE} indicates observed ordinal entries.
#' @param f_sigma_oo_inv_z A function computing The matrix-vector product \eqn{Sigma_{obs, obs}^{-1} * z_{obs}}
#' @param sigma_oo_inv_diag The diagonal of \eqn{Sigma_{obs, obs}^{-1}}
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{Mean for observed ordinal}
#'   \item{\code{var}}{Var for observed ordinal}
#' }
#' @export
#' @keywords internal
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
      "z_new = z
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

        out_trunc = moments_truncnorm(cond_mean_j, cond_std_j, lower[j_in_ord], upper[j_in_ord])
        new_mean_j = out_trunc[['mean']]
        new_var_j = out_trunc[['std']]^2

        if (is.finite(new_mean_j)) z_new[j] = new_mean_j
        if (is.finite(new_var_j)) var_ordinal[j] = new_var_j
      }
      z = z_new"
      new_std = sqrt(1/sigma_oo_inv_diag[ord_in_obs_iter])
      new_mean = z[ord_obs_iter] - (new_std^2) * sigma_oo_inv_z[ord_in_obs_iter]
      a = lower[obs_in_ord_iter]
      b = upper[obs_in_ord_iter]
      out_trunc = moments_truncnorm_vec(mu = new_mean, std = new_std,
                                        a = a, b = b)
      mean_ = out_trunc$mean
      std_ = out_trunc$std
      old_mean = z[ord_obs_iter]
      loc = is.infinite(mean_)
      mean_[loc] = old_mean[loc]
      z[ord_obs_iter] = mean_
      std_[is.infinite(std_)] = 0
      var_ordinal[ord_obs_iter] = std_^2
    }

  }

  list('mean'=z, 'var'=var_ordinal)
}


#' Compute the observed ordinal mean and cov in a row
#'
#' @description  Estimate the conditional mean and cov at observed ordinal entries through sampling
#' @inheritParams update_z_row_ord
#' @param sigma_oo \eqn{Sigma_{obs, obs}}
#' @inheritParams latent_operation_row
#' @param ord_indices Boolean vector where \code{TRUE} indicates ordinal entries.
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{Mean for observed ordinal}
#'   \item{\code{cov}}{Cov for observed ordinal}
#' }
#' @export
#' @keywords internal
est_z_row_ord <- function(z, lower, upper, sigma_oo,
                          ord_indices = NULL, obs_indices = NULL,
                          ord_obs_indices  = NULL, obs_in_ord = NULL, ord_in_obs = NULL,
                          n_sample=5000){
  # need to determine the conditional mean and cov of observed ordinal given observed continuous
  if (is.null(obs_indices)) obs_indices = !is.na(z)
  if (is.null(ord_obs_indices) | is.null(obs_in_ord) | is.null(ord_in_obs)){
    if (is.null(ord_indices)) stop('provide ord_indices')
    ord_obs_indices = ord_indices & obs_indices
    ord_in_obs = ord_obs_indices[obs_indices]
    obs_in_ord = ord_obs_indices[ord_indices]
  }

  p = length(obs_indices)

  if (sum(obs_indices)>1 & any(ord_obs_indices)){
    sigma_ord_ord = sigma_oo[ord_in_obs, ord_in_obs, drop=FALSE]
    cont_in_obs = !ord_in_obs


    # whether to conditional on the continuous observation
    if (!is.null(z) & any(cont_in_obs)){
      z_obs = z[obs_indices]
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
                               n_sample=n_sample)
    if (is.null(z)) z_new = rep(NA, p) else z_new = z
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


#' Compute multivariate truncated normal mean and covariance
#'
#' @description  A wrapper for \code{est_z_row_ord()} where there is no continuous dimension
#' @inheritParams est_z_row_ord
#' @inheritParams latent_operation_row
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{Mean for observed ordinal}
#'   \item{\code{cov}}{Cov for observed ordinal}
#' }
#' @keywords internal
est_z_row_ord_nocont <- function(lower, upper, corr=NULL,
                                 obs_indices = NULL,
                                 cat_input = NULL,
                                 n_sample=5000){
  p = length(lower)
  sigma = corr
  if (is.null(obs_indices)) obs_indices = !is.na(lower)
  z = lower
  if (any(obs_indices)){
    nobs = sum(obs_indices)
    if (is.null(sigma)) sigma_oo = diag(nrow = nobs, ncol = nobs)
    else sigma_oo = sigma[obs_indices,obs_indices,drop=FALSE]

    if (!is.null(cat_input)){
      x_cat = cat_input[['x_cat']]
      cat_index_list = cat_input[['cat_index_list']]
      cat_obs = !is.na(x_cat)
      # A should have dim |cat_o| * |cat_o|
      A = x_to_A(x = x_cat[cat_obs], cat_index_list = cat_index_list[cat_obs], d_cat = sum(obs_indices), old = cat_input[['old']])
      # assuming A satisfying A %*% A = I
      if (!is.null(A)) sigma_oo = A_sigma_tA_at_cat(sigma_oo, A)
    }else A = NULL

    out = est_z_row_ord(z = NULL, lower = lower, upper = upper,
                        sigma_oo = sigma_oo,
                        ord_indices = rep(TRUE, p), obs_indices = obs_indices,
                        n_sample = n_sample)
    mean_est = out$mean[obs_indices]
    if (!is.null(A)) z[obs_indices] = A %*% mean_est else z[obs_indices] = mean_est
  }

  z
}

sample_z_row_ord <- function(z, lower, upper, sigma_oo,
                             ord_indices = NULL, obs_indices = NULL,
                             ord_obs_indices  = NULL, obs_in_ord = NULL, ord_in_obs = NULL,
                             n_sample=5000){
  # need to determine the conditional mean and cov of observed ordinal given observed continuous
  if (is.null(obs_indices)) obs_indices = !is.na(z)
  if (is.null(ord_obs_indices) | is.null(obs_in_ord) | is.null(ord_in_obs)){
    if (is.null(ord_indices)) stop('provide ord_indices')
    ord_obs_indices = ord_indices & obs_indices
    ord_in_obs = ord_obs_indices[obs_indices]
    obs_in_ord = ord_obs_indices[ord_indices]
  }

  p = length(obs_indices)

  if (sum(obs_indices)>1 & any(ord_obs_indices)){
    sigma_ord_ord = sigma_oo[ord_in_obs, ord_in_obs, drop=FALSE]
    cont_in_obs = !ord_in_obs


    # whether to conditional on the continuous observation
    # TODO: only input zcont instead of z
    if (!is.null(z) & any(cont_in_obs)){
      z_obs = z[obs_indices]
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
    "out_ = get_trunc_2dmoments(c(cond_mean), cond_cov,
                               lower_use, upper_use,
                               n_sample=n_sample)"
    zsample = TruncatedNormal::rtmvnorm(n = n_sample,
                                        mu = c(cond_mean), sigma = cond_cov,
                                        lb = lower_use, ub = upper_use)
    # stopifnot(dim(zsample) == c(n_sample, sum(ord_obs_indices)))
    # return from this step for other purposes: sampling or moments estimation
    # record samples
    z_new = matrix(z, n_sample, p, byrow = TRUE)
    z_new[,ord_obs_indices] = zsample
  }else{
    z_new = matrix(z, n_sample, p, byrow = TRUE)
  }

  z_new
}
