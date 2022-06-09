#' Transform a continuous vector to ordinal vector through cutoff function
#'
#' @description  Discretize continuous x to ordinal z with \code{k} levels.
#' @param x continuous vector
#' @param k number of ordinal levels. Generated ordinal vector takes value from \code{1,\ldots,k}
#' @param by Select cut points by absolute values if `dist` and by quantiles if `quantile`
#' @param qmin Cutoff points are slected in the range of `qmin` and `qmax` quantiles of the data `x`
#' @param qmax See `qmin`.
#' @return an ordinal vector.
#' @keywords internal
continuous2ordinal_ = function(x, k = 2, by = 'dist', qmin = 0.1, qmax = 0.9){
  std_dev = sd(x, na.rm = TRUE)
  if (by == 'dist'){
    cutoff = seq(min(x), max(x), length.out = k+1)[-c(1,k+1)]
  }else{
    select = (x>quantile(x,qmin)) & (x<quantile(x,qmax))
    x_sample = x[select]
    if (length(x_sample)<k) stop('Cannot ordinalize the variable')
    cutoff = sample(x_sample, k-1, FALSE)
  }
  cutoff = c(min(x)-0.01, cutoff, max(x)+0.01)
  x = cut(x, cutoff, labels = FALSE, include.lowest = TRUE)
  x
}

#' Transform a continuous vector to ordinal vector through cutoff function.
#'
#' @description  Discretize continuous x to ordinal z with \code{k} levels. If the generated values are guaranteed to have \code{k} unique levels, otherwise, error will be raised.
#' @param x continuous vector
#' @param k number of ordinal levels. Generated ordinal vector takes value from \code{1,\ldots,k}
#' @param by Select cut points by absolute values if `dist` and by quantiles if `quantile`
#' @param qmin Cutoff points are slected in the range of `qmin` and `qmax` quantiles of the data `x`
#' @param qmax See `qmin`.
#' @param max_try The maximum number of attempts to select random cutoffs, since some may yield to ordinal variable with fewer than \code{k} levels
#' @return an ordinal vector.
#' @export
#' @keywords internal
continuous2ordinal <- function(x, k = 2, by = 'dist', qmin = 0.1, qmax = 0.9, max_try=20){
  success = FALSE
  for (i in 1:max_try){
    result = continuous2ordinal_(x,k,by=by,qmin=qmin,qmax=qmax)
    if (length(unique(result))==k){
      success = TRUE
      break
    }
  }
  result
}

#' Sample random correlation matrix
#'
#' @description Sample full rank random correlation matrix
#' @param p Data dimension
#' @return a \code{p*p} correlation matrix
#' @export
generate_sigma <- function(p){
  W = matrix(rnorm(p*p), nrow = p)
  Sigma = W %*% t(W)
  cov2cor(Sigma)
}

#' Wrapper
#' @export
#' @keywords internal
generate_mixed_from_gc <- function(...) generate_gc(...)

#' Sample data from Gaussian copula distribution
#'
#' @description  Sample mixed data vector containing continuous, ordinal and binary data from a Gaussian copula model
#' @param rank Only used for \code{generate_lrgc}
#' @param n number of samples to draw
#' @param seed random seed
#' @param var_types A list indicating the locations of continuous, ordinal and binary variables. The sum of its elements, \code{p}, is the desired number of variables
#' @param data_only If \code{TRUE}, only return the generated data matrix
#' @param cont_transform A monotonic function to be applied to continuous columns
#' @param cutoff_by When ordinalizing, select cut points by absolute values if \code{dist} and by quantiles if \code{quantile}
#' @param qmin Cutoff points are selected in the range of \code{qmin} and \code{qmax} quantiles of the data
#' @param qmax See \code{qmin}
#' @param num_ord Number of ordinal levels to use
#' @param corr Only used for \code{generate_gc}. \code{NULL}, or a correlation matrix or a list of copula correlation matrices
#' @param W Only used for \code{generate_lrgc}.  \code{NULL}, or a matrix of dimension \code{p} by \code{rank}
#' @param sigma Only used for \code{generate_lrgc}. A scalar between 0 and 1.
#' @details For \code{generategc}, if \code{corr} is \code{NULL}, it will be allocated a random correlation matrix. For \code{generate_lrgc},  if \code{W} is \code{NULL}, it will be allocated a random matrix. \code{corr} used in LRGC is \code{Wt(W)+sigma*I}. If \code{corr} is a single correlation matrix, \code{n} samples will be generated from the Gaussian copula model with \code{corr} as the copula correlation. If \code{corr} is a list of correlation matrices, for each element of \code{corr}, \code{n} samples will be generated from the Gaussian copula model with that element as its copula correlation.
#' @name generate_copula
NULL
#> NULL

#' @describeIn generate_copula Generate samples from a full rank Gaussian copula
#' @export
#' @examples
#' generate_gc(n=500, var_types = list('cont'=1:5, 'ord'=6:10, 'bin'=11:15))
#' generate_gc(n=500, var_types = list('cont'=1:5, 'ord'=6:10), cont_transform=function(x)x^3)
#' generate_gc(n=500, var_types = list('cont'=1:5, 'ord'=6:10), num_ord=3)
generate_gc <- function(n=2000,  var_types = NULL,
                        cont_transform=NULL,
                        cutoff_by='quantile', qmin=0.05, qmax=0.95, num_ord=5,
                        corr=NULL, seed=NULL, data_only=TRUE){
  sigma = corr
  if (is.null(cont_transform)) cont_transform <- function(x) qexp(pnorm(x), rate = 1/3)
  if (is.null(var_types)) var_types = list('cont'=1:5, 'ord'=6:10, 'bin'=11:15)
  cont_index = var_types$cont
  ord_index = var_types$ord
  bin_index = var_types$bin
  all_index = c(cont_index, ord_index, bin_index)
  p = length(all_index)
  if (min(all_index)!= 1 | max(all_index) != p | length(unique(all_index))!=p){
    stop('Inconcistent specification of variable types indexing')
  }

  if (!is.null(seed)) set.seed(seed)
  if (is.null(sigma)) sigma = generate_sigma(p)
  if (!is.list(sigma)) sigma = list(sigma)
  if (ncol(sigma[[1]])!=p) stop('Inconcistent dimension between variable lengths and copula correlation')

  l = length(sigma)
  X = vector('list', l)
  for (i in 1:l){
    X[[i]] <- MASS::mvrnorm(n, mu = numeric(p), Sigma = sigma[[i]])
  }
  X = do.call(rbind, X)

  X[,cont_index] = cont_transform(X[,cont_index])
  for (i in ord_index) X[,i] = continuous2ordinal(X[,i], k=num_ord, by = cutoff_by, qmin = qmin, qmax = qmax)
  for (i in bin_index) X[,i] = continuous2ordinal(X[,i], k=2, by = cutoff_by, qmin = qmin, qmax = qmax)

  if (data_only) X else list(X = X, corr = sigma)
}

#' @describeIn generate_copula Generate samples from a low rank Gaussian copula
#' @export
#' @examples
#' generate_lrgc(rank = 5, sigma = 0.1, n=500, var_types = list('cont'=1:100), cont_transform=function(x)x^3)
#' generate_lrgc(rank = 5, sigma = 0.1, n=500, var_types = list('ord'=1:100), num_ord=10)
#' generate_lrgc(rank = 5, sigma = 0.1, n=500, var_types = list('bin'=1:100))
generate_lrgc <- function(rank=5, sigma = 0.1, n=2000, var_types = NULL,
                          cont_transform=NULL, cutoff_by='quantile', qmin=0.05, qmax=0.95, num_ord=5,
                          W=NULL,seed=NULL,data_only=TRUE){
  if (is.null(cont_transform)) cont_transform <- function(x) qexp(pnorm(x), rate = 1/3)
  if (is.null(var_types)) var_types = list('cont'=1:20, 'ord'=21:40, 'bin'=41:60)
  cont_index = var_types$cont
  ord_index = var_types$ord
  bin_index = var_types$bin
  all_index = c(cont_index, ord_index, bin_index)
  p = length(all_index)
  if (min(all_index)!= 1 | max(all_index) != p | length(unique(all_index))!=p){
    stop('Inconcistent specification of variable types indexing')
  }

  if (!is.null(seed)) set.seed(seed)
  if (is.null(W)) W = matrix(rnorm(p*rank), nrow = p)
  else stopifnot(dim(W) != c(p,rank))
  stopifnot(sigma>0 & sigma<1)
  if (!is.list(W)) W = list(W)


  l = length(W)
  X = vector('list', l)
  for (i in 1:l){
    W[[i]] = t(apply(W[[i]], 1, function(x){x/sqrt(sum(x^2))  * sqrt((1-sigma))}))
    X[[i]] <- matrix(rnorm(n*rank), ncol = rank) %*% t(W[[i]]) + matrix(rnorm(n*p, sd = sqrt(sigma)), ncol = p)
  }
  X = do.call(rbind, X)

  X[,cont_index] = cont_transform(X[,cont_index])
  for (i in ord_index) X[,i] = continuous2ordinal(X[,i], k=num_ord, by = cutoff_by, qmin = qmin, qmax = qmax)
  for (i in bin_index) X[,i] = continuous2ordinal(X[,i], k=2, by = cutoff_by, qmin = qmin, qmax = qmax)

  if (data_only) X else list(X = X, W = W, sigma = sigma)
}

gen_nominal_copula_corr <- function(p_cat_vec, p_noncat, seed=NULL, eps=1e-2){
  if (!is.null(seed)) set.seed(seed)
  last = 0
  l =  length(p_cat_vec)
  cat_index_list = vector('list', l)
  for (i in 1:l){
    cat_index_list[[i]] = last + (1:p_cat_vec[i])
    last = last + length(cat_index_list[[i]])
  }
  p_catall = sum(p_cat_vec)
  p = p_catall + p_noncat
  sigma = gcimputeR::generate_sigma(p = p)
  #
  "  A = diag(p)
  for (i in 1:l){
    index = cat_index_list[[i]]
    o_eigen = eigen(sigma[index,index])
    m = diag(1/sqrt(o_eigen$values)) %*% t(o_eigen$vectors)
    A[index,index] = m
  }
  sigma = A %*% sigma %*% t(A)
  for (i in 1:l){
    index = cat_index_list[[i]]
    sigma[index,index] = diag(length(index))
  }"
  sigma = project_to_nominal_corr(sigma, cat_index_list)
  if (min(eigen(sigma)$values)<eps){
    sigma = cov2cor(sigma + diag(p) * eps)
  }
  sigma
}

gen_nominal_copula <- function(n, p_cat_vec, p_noncat, p_ord=0,
                               mu_scale = 0.5, seed=NULL,
                               mu_cat = NULL, sigma=NULL, trans_cont=NULL,
                               ord_num = 5,
                               old = FALSE){
  if (old) return(gen_nominal_copula_old(n, p_cat_vec, p_noncat, mask_fraction=0,
                                         mu_scale, seed, mu_cat, sigma))
  #print(mu_scale)
  if (is.null(trans_cont)) trans_cont = function(x) x
  if (!is.null(seed)) set.seed(seed)
  l = length(p_cat_vec)
  cat_index_list = vector('list', l)
  last = 0
  starts = numeric(l)
  for (i in 1:l){
    cat_index_list[[i]] = last + (1:p_cat_vec[i])
    starts[i] = last + 1
    last = last + length(cat_index_list[[i]])
  }
  p_cat = sum(p_cat_vec)
  l_cat = length(p_cat_vec)
  p = p_cat + p_noncat
  #
  if (is.null(mu_cat)){
    mu_cat = rnorm(p_cat)*mu_scale
    mu_cat[starts] = 0
  }
  if (is.null(sigma)) sigma = gen_nominal_copula_corr(p_cat_vec, p_noncat)
  z = MASS::mvrnorm(n = n, mu = c(mu_cat, rep(0, p_noncat)), Sigma = sigma)
  x = matrix(0, n, l_cat+p_noncat)
  for (i in 1:l){
    index = cat_index_list[[i]]
    x[,i] = apply(z[,index], 1, which.max)
  }
  if (p_ord>0){
    for (j in 1:p_ord) x[,l_cat+j] = continuous2ordinal(z[,p_cat+j], k=ord_num)
    if (p_ord < p_noncat) for (j in (p_ord+1):p_noncat) x[,l_cat+j] = trans_cont(z[,p_cat+j])
  }else{
    x[,-(1:l_cat)] = trans_cont(z[,-(1:p_cat)])
  }
  #if (mask_fraction>0) xmask = mask_MCAR(x, mask_fraction) else xmask = x
  #df = as.data.frame(xmask)
  #for (i in 1:l) df[,i] = as.factor(df[,i])
  #xmask_cat = xmask[,1:l,drop=FALSE]
  #zmask_noncat = xmask[,-(1:l),drop=FALSE]
  # 'xmask'=xmask, 'df_xmask'=df,'xmask_cat'=xmask_cat,'zmask_noncat'=zmask_noncat,
  list('x'=x,
       'sigma'=sigma, 'mu_cat'=mu_cat, 'z'=z,
       'cat_index_list'=cat_index_list,
       'cat_index'=1:l_cat
  )
}

gen_nominal_copula_old <- function(n, p_cat_vec, p_noncat, mask_fraction=0, mu_scale = 0.5, seed=NULL,
                                   mu_cat = NULL, sigma=NULL){
  if (!is.null(seed)) set.seed(seed)
  l = length(p_cat_vec)
  cat_index_list = vector('list', l)
  last = 0
  for (i in 1:l){
    cat_index_list[[i]] = last + (1:p_cat_vec[i])
    last = last + length(cat_index_list[[i]])
  }
  p_cat = sum(p_cat_vec)
  l_cat = length(p_cat_vec)
  p = p_cat + p_noncat
  #
  if (is.null(mu_cat)) mu_cat = rnorm(p_cat)*0.5
  if (is.null(sigma)) sigma = gen_nominal_copula_corr(p_cat_vec, p_noncat)
  z = MASS::mvrnorm(n = n, mu = c(mu_cat, rep(0, p_noncat)), Sigma = sigma)
  x = matrix(0, n, l_cat+p_noncat)
  for (i in 1:l){
    index = cat_index_list[[i]]
    zuse = cbind(matrix(0,n,1), z[,index])
    x[,i] = apply(zuse, 1, which.max)
  }
  x[,-(1:l_cat)] = z[,-(1:p_cat)]
  if (mask_fraction>0) xmask = mask_MCAR(x, mask_fraction) else xmask = x
  df = as.data.frame(xmask)
  for (i in 1:l) df[,i] = as.factor(df[,i])

  xmask_cat = xmask[,1:l,drop=FALSE]
  zmask_noncat = xmask[,-(1:l),drop=FALSE]
  list('x'=x, 'xmask'=xmask, 'df_xmask'=df,
       'xmask_cat'=xmask_cat,
       'zmask_noncat'=zmask_noncat,
       'sigma'=sigma, 'mu_cat'=mu_cat, 'z'=z,
       'cat_index_list'=cat_index_list,
       'cat_index'=1:l_cat
  )
}

