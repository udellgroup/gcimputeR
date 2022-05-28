#' Transform original observation X to copula observation Z
#'
#' @description  For continuous \code{X}, compute its corresonding latent value \code{Z}; For ordinal \code{X}, compute its corresponding truncated interval on latent value \code{Z}.
#' @param X Original observation matrix
#' @param type Data type. Either \code{"continuous"} or \code{"ordinal"}.
#' @return A list containing
#' \describe{
#'   \item{\code{Z}}{Available when \code{type=continuous}. Corresponding value}
#'   \item{\code{Lower}}{Available when \code{type=ordinal}. Lower boundary of corresponding truncated interval}
#'   \item{\code{Upper}}{Available when \code{type=ordinal}. Upper boundary of corresponding truncated interval}
#' }
#' @export
range_transform = function(X, type){
  n = dim(X)[1]
  p = dim(X)[2]

  # for continuous dimensions, return a transformed value for each observation
  if (type == 'continuous' & p>0){
    r_val = X
    for (j in 1:p){
      fun = ecdf(X[,j])
      r_val[,j] = fun(X[,j])
    }
    r_val = qnorm(r_val * n/(n+1))
    return(list(Z = r_val))
  }

  # for ordinal dimensions, return a transformed interval for each observation
  if (type == 'ordinal' & p>0){
    d = numeric(p)
    r_lower = X
    r_upper = X
    for (j in 1:p){
      x = sort(unique(X[,j]))
      if (length(x) <= 1) stop(paste('only one level in dimension', j))
      d[j] = min(x[-1] - x[-length(x)])/2

      fun = ecdf(X[,j])
      r_lower[,j] = fun(X[,j] - d[j])
      r_upper[,j] = fun(X[,j] + d[j])
    }
    #r_lower = qnorm(r_lower * n/(n+1))
    #r_upper = qnorm(r_upper * n/(n+1))
    r_lower = qnorm(r_lower)
    r_upper = qnorm(r_upper)
    return(list(Lower=r_lower, Upper=r_upper))
  }
}



#' Transform imputed Z matrix to X in the observed space
#'
#' @description  Impute orginal incomplete observation using imputed \code{Z} matrix and marginal distribution information drawn from original observation \code{X}
#' @param Z imputed Z matrix
#' @param X original observation matrix
#' @param d_index ordinal dimension indexes.
#' @return Imputed data matrix
#' @export
Ximp_transform = function(Z, X, d_index){
  # for continuous, use default empirical quantile
  n = dim(Z)[1]
  p = dim(Z)[2]
  Ximp = X
  c_index = !d_index

  # ordinal
  for (j in which(d_index)){
    miss_ind = is.na(X[,j])
    n_obs = n-sum(miss_ind)
    # decimal to integer 1 to n_obs
    xmis_loc = pmax(ceiling(pnorm(Z[miss_ind,j]) * n_obs), 1)
    # If extremely small value appear in Z[miss_ind, ] such that xmis_loc contains 0,
    # replace it with 1 (impute using the smallest observation)
    Ximp[miss_ind,j] = sort(X[!miss_ind,j])[xmis_loc]
  }

  # continuous
  for (j in which(c_index)){
    miss_ind = is.na(X[,j])
    Ximp[miss_ind,j] = quantile(X[!miss_ind,j], pnorm(Z[miss_ind,j]))
  }

  Ximp
}

Ximp_transform_cont = function(Z, X){
  # for continuous, use default empirical quantile
  n = dim(Z)[1]
  p = dim(Z)[2]
  stopifnot(p>=1)
  Ximp = X

  for (j in 1:p){
    miss_ind = is.na(X[,j])
    Ximp[miss_ind,j] = quantile(X[!miss_ind,j], pnorm(Z[miss_ind,j]))
  }

  Ximp
}

Ximp_transform_ord = function(Z,X){
  n = dim(Z)[1]
  p = dim(Z)[2]
  stopifnot(p>=1)
  Ximp = X

  # ordinal
  for (j in 1:p){
    miss_ind = is.na(X[,j])
    n_obs = n-sum(miss_ind)
    # decimal to integer 1 to n_obs
    xmis_loc = pmax(ceiling(pnorm(Z[miss_ind,j]) * n_obs), 1)
    # If extremely small value appear in Z[miss_ind, ] such that xmis_loc contains 0,
    # replace it with 1 (impute using the smallest observation)
    Ximp[miss_ind,j] = sort(X[!miss_ind,j])[xmis_loc]
  }

  Ximp
}

Ximp_transform_cat <- function(Z_cat, X_cat, cat_index_list, old = FALSE){
  if (ncol(Z_cat) != sum(purrr::map_int(cat_index_list, length))) stop('something wrong')
  cat_index_list = adjust_index_list(cat_index_list)
  Ximp_cat = X_cat
  for (j in seq_along(cat_index_list)){
    index_m = is.na(X_cat[,j])
    index_cat = cat_index_list[[j]]
    zmis = Z_cat[index_m,index_cat,drop=FALSE]
    Ximp_cat[index_m,j] = apply(zmis, 1, nominal_z_to_x_col, old = old)
  }
  Ximp_cat
}

#' From observed X to latent Z
#'
#' @description  Given data matrix \code{X}, prepare the latent matrix \code{Z} and the value bounds for the ordinal dimensions
#' @inheritParams range_transform
#' @inheritParams initZ_interval_truncated
#' @param d_index Boolean vector with \code{TRUE} at ordinal dimensions
#' @inheritParams initZ_noncat
#' @return A list containing
#' \describe{
#'   \item{\code{Z}}{Transformed latent matrix}
#'   \item{\code{Lower}}{Lower boundary for ordinal dimensions. \code{NA} at missing locations.}
#'   \item{\code{Upper}}{Upper boundary for ordinal dimensions. \code{NA} at missing locations.}
#' }
observed_to_latent <- function(X, d_index, method='univariate_mean'){
  c_index = !d_index
  Z = matrix(NA, nrow(X), ncol(X))
  if (any(d_index)){
    r = range_transform(X[,d_index,drop=FALSE], type= 'ordinal')
    Z_ord_lower = r$Lower
    Z_ord_upper = r$Upper
    Z_ord = initZ_interval_truncated(Z_ord_lower, Z_ord_upper, method=method)
    Z[, d_index] = Z_ord
  }else{
    Z_ord_lower = NULL
    Z_ord_upper = NULL
  }

  if (any(c_index)){
    Z_cont = range_transform(X[,c_index,drop=FALSE], type = 'continuous')$Z
    Z[, c_index] = Z_cont
  }else{
    Z_cont = NULL
  }

  list(Z=Z, Lower=Z_ord_lower, Upper=Z_ord_upper)
}

latent_to_observed <- function(Zimp, X, mu, cat_labels, ord_in_noncat,
                               cat_index, cat_index_all, cat_index_list, old=FALSE){
  d_cat = length(cat_index_all)
  n = nrow(Zimp)
  Z_cat = Zimp[,cat_index_all] + matrix(mu, n, d_cat, byrow = TRUE)
  Ximp = X
  X_cat = X[,cat_index,drop=FALSE]
  X_noncat = X[,!cat_index,drop=FALSE]
  Ximp[,cat_index] = Ximp_transform_cat(Z_cat = Z_cat, X_cat = X_cat,
                                        cat_index_list = cat_index_list, old = old) #***!
  Ximp[,!cat_index] = Ximp_transform(Z = Zimp[,-cat_index_all], X = X_noncat, d_index = ord_in_noncat)

  cat_index_int = which(cat_index)
  for (j in cat_index_int) Ximp[,j] = relabel(Ximp[,j], cat_labels[[as.character(j)]])

  Ximp
}

#' Get the lower and upper bounds for truncated normal moments calculation at categorical columns
#'
#' @description  Impute orginal incomplete observation using imputed \code{Z} matrix and marginal distribution information drawn from original observation \code{X}
#' @param X_cat (incomplete) categorical data matrix
#' @param mu mean vector for categorical
#' @param cat_index_list list indicating the number of levels for each categorical variable
#' @return Imputed data matrix
get_cat_bounds <- function(X_cat, mu, cat_index_list, check=FALSE, old = FALSE){
  # TODO check values of X_cat
  if (old) return(get_cat_bounds_old(X_cat, mu, cat_index_list, check))
  d_cat = sum(purrr::map_int(cat_index_list, length))
  if (d_cat != length(mu)) stop('invalid input')
  p_cat = length(cat_index_list)
  if (p_cat != ncol(X_cat)) stop('invalid input')
  n = nrow(X_cat)
  lower = matrix(NA, n, d_cat)
  upper = matrix(NA, n, d_cat)
  incat_index_list = adjust_index_list(cat_index_list)
  for (j in 1:p_cat){
    index_o = !is.na(X_cat[,j])
    x_cat_obs = X_cat[index_o,j]
    n_obs = length(x_cat_obs)
    # initialize to (0, Inf): at argmax, we want (-inf, inf), at other loc, we want (mu_j - mu_argmax, inf)
    # index in 1,...,d_cat
    index_cat = incat_index_list[[j]]
    dj_cat = length(index_cat)
    l_o = matrix(0, n_obs, dj_cat)
    u_o = l_o + Inf

    # adjust for mean
    # length d_cat
    mu_j = mu[index_cat]
    # z_argmax - z_{-argmax} + mu_j[argmax] - mu_j[-argmax] >=0
    # thus z_argmax - z_{-argmax} >= mu_j[-argmax] - mu_j[argmax] (RHS computed below)
    mu_diff = matrix(mu_j, n_obs, dj_cat, byrow = TRUE)  - matrix(mu_j[x_cat_obs], n_obs, dj_cat, byrow = FALSE)
    l_o = l_o + mu_diff
    # no constraints at argmax, thus -Inf lower
    argmax_coor = matrix(c(1:n_obs, x_cat_obs), nrow = n_obs)
    l_o[argmax_coor] = -Inf

    lower[index_o, index_cat] = l_o
    upper[index_o, index_cat] = u_o
  }

  if (check){
    for (i in 1:n){
      ind1 = is.na(lower[i,])
      ind2 = get_cat_slicing_index(X_cat[i,], incat_index_list, keep = 'missing', d_cat=d_cat)$cat
      if (!all(ind1==ind2)) stop('something wrong!')
    }
  }
  list(lower = lower, upper = upper)
}

get_cat_bounds_old <- function(X_cat, mu, cat_index_list, check=FALSE){
  d_cat = sum(purrr::map_int(cat_index_list, length))
  if (d_cat != length(mu)) stop('invalid input')
  p_cat = length(cat_index_list)
  if (p_cat != ncol(X_cat)) stop('invalid input')
  n = nrow(X_cat)
  lower = matrix(NA, n, d_cat)
  upper = matrix(NA, n, d_cat)
  incat_index_list = adjust_index_list(cat_index_list)
  for (j in 1:p_cat){
    index_o = !is.na(X_cat[,j])
    x_cat_obs = X_cat[index_o,j]
    n_obs = length(x_cat_obs)
    # adjust for transforming hyperplane transformation into hypercube transformation
    # bounds computed after this step is a special case when mu = 0
    index_cat = incat_index_list[[j]]
    dj_cat = length(index_cat)
    l_o = matrix(0, n_obs, dj_cat)
    u_o = l_o + Inf
    index_base = x_cat_obs == 1
    l_o[index_base,] = -Inf
    u_o[index_base,] = 0

    # adjust for mean
    mu_use = matrix(mu[index_cat], n_obs, dj_cat, byrow = TRUE)
    row_index = which(!index_base)
    col_index = x_cat_obs[row_index]-1
    select = matrix(c(row_index, col_index), nrow = length(row_index))
    mu_use[!index_base,] = mu_use[select] - mu_use[!index_base,]
    mu_use[select] = mu[index_cat][col_index]

    lower[index_o, index_cat] = l_o - mu_use
    upper[index_o, index_cat] = u_o - mu_use
  }

  if (check){
    for (i in 1:n){
      ind1 = is.na(lower[i,])
      ind2 = get_cat_slicing_index(X_cat[i,], cat_index_list, keep = 'missing', d_cat=d_cat)$cat
      if (!all(ind1==ind2)) stop('something wrong!')
    }
  }
  list(lower = lower, upper = upper)
}

