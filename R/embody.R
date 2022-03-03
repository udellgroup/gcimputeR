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

est_z_row_ord <- function(z_obs, lower, upper, obs_indices, ord_obs_indices, ord_in_obs, sigma_oo,
                          method='Explicit', n_sample=5000){
  # need to determine the conditional mean and cov of observed ordinal given observed continuous
  if (sum(obs_indices)>1 & any(ord_obs_indices)){
    cont_in_obs = !ord_in_obs
    sigma_ord_ord = sigma_oo[ord_in_obs, ord_in_obs]
    sigma_cont_ord = sigma_oo[cont_in_obs, ord_in_obs]
    sigma_cont_cont = sigma_oo[cont_in_obs, cont_in_obs]

    tot_m = cbind(matrix(z_obs[cont_in_obs], ncol = 1), sigma_cont_ord)
    sol_m = solve(sigma_cont_cont, tot_m)

    cond_mean = t(sigma_cont_ord) %*% sol_m[,1]
    cond_cov = sigma_ord_ord - t(sigma_cont_ord) %*% sol_m[,-1]

    out = get_trunc_2dmoments(cond_mean, cond_cov, lower, upper, method=method, n_sample=n_sample)
  }else{
    out = NULL
  }

  # return: list with names 'mean' and 'cov' or NULL
  out
}

get_trunc_2dmoments <- function(mean, cov, lower, upper, method='Explicit', n_sample=5000){
  if (length(mean) != ncol(cov)){
    print(length(mean))
    print(dim(cov))
    stop('inconsistent shapes')
  }
  p = length(mean)
  switch (method,
          'TruncatedNormal' = {
            z = TruncatedNormal::rtmvnorm(n = n_sample, mu = mean, sigma = cov, lb = lower, ub = upper)
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
