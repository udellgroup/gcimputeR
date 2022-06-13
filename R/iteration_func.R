#' Operation to conduct in the latent space
#'
#' @description  Conduct one of three tasks in the latent space. If \code{task = "em"}, conduct an em iteration. If \code{task = "fillup"}, impute the missing entries of \code{Z} as their conditional mean. If \code{task = "sample"}, conduct multiple imputation on the missing entries of \code{Z}.
#' @param task Task to perform. One of \code{"em", "fillup", "sample"}.
#' @param corr Current copula correlation estimate
#' @inheritParams em_mixedgc
#' @inheritParams impute_mixedgc
#' @return A list containing
#' \describe{
#'   \item{\code{corr}}{Available when \code{task = 'em'}. Updated correlation estimate.}
#'   \item{\code{loglik}}{Available when \code{task = 'em'}. The average log-likelihood.}
#'   \item{\code{Z}}{Available when \code{task = 'em'}. Incomplete \code{Z} with updated observed ordinal entries}
#'   \item{\code{Zimp}}{Available when \code{task = 'em'} or \code{task == 'fillup'} . Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
#'   \item{\code{Zimp_sample}}{Available when \code{task = 'sample'}. Multiple imputation samples.}
#'   \item{\code{C}}{Available when \code{task = 'em'}. The conditional co-variance due to missingness}
#'   \item{\code{var_ordinal}}{Available when \code{task = 'em'} or \code{task = 'fillup'}. The conditional variance due to truncation, i.e. Var(z|a < z < b)}
#' }
#' @export
#' @keywords internal
latent_operation <- function(task,
                             Z, Lower, Upper,
                             d_index, dcat_index,
                             cat_input,
                             corr,
                             trunc_method='Iterative', n_update=1, n_sample=5000,
                             n_MI = 1){
  n = nrow(Z)
  p = ncol(Z)
  # TODO: remove var alias
  Z_lower = Lower
  Z_upper = Upper
  if (!is.logical(d_index)) stop('must use logical d_index')

  switch (task,
          'em' = {
            out = list('Z'=Z, 'Zimp'=Z, 'loglik'=0, 'C'=matrix(0,p,p), 'var_ordinal'=matrix(0,n,p))
          },
          'fillup' = {
            out = list('Zimp'=Z,  'var_ordinal'=matrix(0,n,p))
          },
          'sample' = {
            out = list('Zimp_sample'= array(0, dim = c(n,p,n_MI)))
          }
  )

  if (!is.null(cat_input)){
    cat_input_row = list(x_cat=NULL,
                         cat_index_list=cat_input$cat_index_list,
                         cat_index_all=cat_input$cat_index_all,
                         old = cat_input$old) #***
  }else cat_input_row = NULL

  for (i in 1:n){
    if (!is.null(cat_input_row)) cat_input_row[['x_cat']] = cat_input$X_cat[i,]
    row_out = latent_operation_row(task,
                                   Z[i,], Z_lower[i,], Z_upper[i,],
                                   d_index=d_index, dcat_index=dcat_index,
                                   cat_input=cat_input_row,
                                   corr=corr,
                                   trunc_method = trunc_method,n_update=n_update, n_sample=n_sample,
                                   n_MI = n_MI)

    switch (task,
            'em' = {
              for (name in c('Z', 'Zimp', 'var_ordinal')){
                out[[name]][i,] = row_out[[name]]
              }
              for (name in c('loglik', 'C')){
                out[[name]] = out[[name]] + row_out[[name]]/n
              }
            },
            'fillup' = {
              for (name in c('Zimp', 'var_ordinal')){
                out[[name]][i,] = row_out[[name]]
              }
            },
            'sample' = {
              out[['Zimp_sample']][i,,] = row_out[['Zimp_sample']]
            }
    )
  }

  if (task == 'em'){
    out[['corr']] = cov(out[['Zimp']]) + out[['C']]
    if (any(is.na(out[['corr']]))) stop('invalid correlation')
  }

  out
}

latent_operation_LRGC <- function(task,
                                  Z, Lower, Upper,
                                  d_index,
                                  U, d, sigma,
                                  dcat_index=NULL,cat_input=NULL,
                                  trunc_method='Iterative', n_update=1, n_sample=5000,
                                  n_MI = 1){
  n = nrow(Z)
  p = ncol(Z)
  rank = ncol(U)
  if (!is.logical(d_index)) stop('must use logical d_index')

  switch (task,
          'em' = {
            out = list('Z'=Z, 'var_ordinal'=matrix(0,n,p), 'loglik'=0,
                       'A'=array(0, dim=c(n,rank,rank)), 'S'=matrix(0,n,rank), 'SS'=array(0, dim=c(n,rank,rank)),
                        'zobs_norm' = 0)
          },
          'fillup' = {
            out = list('Zimp'=Z,  'var_ordinal'=matrix(0,n,p))
          },
          'sample' = {
            out = list('Zimp_sample'= array(0, dim = c(n,p,n_MI)))
          }
  )

  if (!is.null(cat_input)){
    cat_input_row = list(x_cat=NULL,
                         cat_index_list=cat_input$cat_index_list,
                         cat_index_all=cat_input$cat_index_all,
                         old = cat_input$old) #***
  }else cat_input_row = NULL

  for (i in 1:n){
    if (!is.null(cat_input_row)) cat_input_row[['x_cat']] = cat_input$X_cat[i,]
    row_out = latent_operation_LRGC_row(task,
                                        Z[i,], Lower[i,], Upper[i,],
                                        d_index=d_index,
                                        U=U, d=d, sigma=sigma,
                                        cat_input=cat_input_row, dcat_index=dcat_index,
                                        trunc_method = trunc_method,n_update=n_update, n_sample=n_sample,
                                        n_MI = n_MI)

    switch (task,
            'em' = {
              for (name in c('Z', 'var_ordinal','S')){
                out[[name]][i,] = row_out[[name]]
              }
              for (name in c('A', 'SS')){
                out[[name]][i,,] = row_out[[name]]
              }
              for (name in c('loglik', 'zobs_norm')){
                out[[name]] = out[[name]] + row_out[[name]]
              }
            },
            'fillup' = {
              for (name in c('Zimp', 'var_ordinal')){
                out[[name]][i,] = row_out[[name]]
              }
            },
            'sample' = {
              out[['Zimp_sample']][i,,] = row_out[['Zimp_sample']]
            }
    )
  }

  out
}

Mstep_LRGC <- function(Z, Cord, U, sigma, A, S, SS){
  W = array(0, dim = dim(U))
  p = ncol(Z)
  s = 0
  for (j in 1:p){
    index_j = !is.na(Z[,j])
    r = Mstep_LRGC_col(Z[index_j,j], Cord[index_j,j],
                       U = U[j,], sigma = sigma,
                       A = A[index_j,,,drop=FALSE],
                       S = S[index_j,,drop=FALSE],
                       SS = SS[index_j,,,drop=FALSE])
    W[j,] = r$w
    s = s + r$s
  }
  return(list(W = W, s = s))
}

Mstep_LRGC_col <- function(Z, Cord, U, sigma, A, S, SS){
  rj = sum_2d_scale(M = S, v = Z)
  rj = rj + sum_3d_scale(A, v = Cord) %*% matrix(U, ncol=1)
  n = dim(A)[1]
  Fj = sum_3d_scale(SS+sigma*A, v = rep(1,n))
  w = solve(Fj,rj)
  s = sum(rj * w)
  return(list(w = w, s = s))
}



#' Each iteration of EM algorithm fit (LRGC)
#'
#' @description  fit the Gaussian copula model with incomplete mixed type data
#' @param Z Incomplete Z matrix
#' @param Z_lower Lower boundary of truncated intervals for ordinal columns
#' @param Z_upper Upper boundary of truncated intervals for ordinal columns
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
#' @keywords internal
em_mixedgc_ppca_iter = function(Z, Lower, Upper, d_index, W, sigma, n_update=1){
  Z_lower = Lower
  Z_upper = Upper
  n = nrow(Z)
  p = ncol(Z)
  rank = ncol(W)

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
    ord_indices = d_index
    mis_indices = is.na(Z[i,])
    obs_indices = !mis_indices
    ord_obs_indices = ord_indices & obs_indices
    ord_in_obs = ord_obs_indices[obs_indices]
    obs_in_ord = ord_obs_indices[ord_indices]

    #d_in_o = which(obs_indices <= k)
    #index_d = obs_indices[d_in_o]

    #
    zi_obs = matrix(Z[i,obs_indices], ncol = 1)
    Ui_obs = matrix(U[obs_indices,], ncol = rank)
    UU_obs = t(Ui_obs) %*% Ui_obs

    # used in both ordinal and factor block
    ans = solve(sigma * diag(d^{-2}) + UU_obs, cbind(diag(rank), t(Ui_obs)))
    AU = ans[,-(1:rank),drop=FALSE]
    Ai = ans[,1:rank,drop=FALSE]
    A[i,,] = Ai


    # Only implement when there is at least one observed ordinal dimension(to be imputed)
    # and another observed dimension (ordinal or continuous) used as information to impute
    f_sigma_oo_inv_z = function(zobs) (zobs - (Ui_obs %*% (AU %*% zobs)))/sigma
    sigma_oo_inv_diag = (1 - quad_mul(Ai, Ui_obs))/sigma
    # truncated moments
    out <- update_z_row_ord(Z[i,], Lower[i,], Upper[i,],
                            obs_indices = obs_indices,
                            ord_obs_indices = ord_obs_indices,
                            ord_in_obs = ord_in_obs,
                            obs_in_ord = obs_in_ord,
                            f_sigma_oo_inv_z = f_sigma_oo_inv_z,
                            sigma_oo_inv_diag = sigma_oo_inv_diag,
                            n_update = n_update)
    Z[i,] = out$mean
    C[i,] = out$var

    zi_obs = matrix(Z[i,obs_indices], ncol = 1)
    si = matrix(AU %*% zi_obs, ncol = 1)
    S[i,] = si
    tSC[i,,] = si %*% t(si)
    SC[i,,] = AU %*% (C[i,obs_indices] * t(AU))

    # pseudo-likelihood
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
    AC[j,,] = sum_3d_scale(A, v = C[,j], index = index_j)
    rj = sum_2d_scale(M = S, v = Z[,j], index = index_j) + AC[j,,] %*% matrix(U[j,], ncol=1)
    # denominator
    ASC[j,,] = sum_3d_scale(SC, v = rep(1,n), index = index_j) + sum_3d_scale(A, v = rep(sigma,n), index = index_j)
    Fj = sum_3d_scale(tSC, v = rep(1,n), index = index_j) +  ASC[j,,]
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

  return(list(W = est$W, sigma = est$sigma, loglik=-negloglik/2, Z=Z, S = S, C=C))
}



