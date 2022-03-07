#' Transform a continuous vector to ordinal vector through cutoff function
#'
#' @description  Discretize continuous x to ordinal z with \code{k} levels.
#' @param x continuous vector
#' @param k number of ordinal levels. Generated ordinal vector takes value from \code{1,\ldots,k}
#' @param by Select cut points by absolute values if `dist` and by quantiles if `quantile`
#' @param qmin Cutoff points are slected in the range of `qmin` and `qmax` quantiles of the data `x`
#' @param qmax See `qmin`.
#' @return an ordinal vector.
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


#' Sample data from Gaussian copula distribution
#'
#' @description  Sample mixed data vector consisting of continuous, ordinal and binary marginals from a Gaussian copula model
#' @param sigma a correlation matrix or a list of copula correlation matrices
#' @param n number of samples to draw for each copula correlaiton matrix
#' @param seed random seed
#' @param var_types A list indicating the locations of continuous, ordinal and binary variables
#' @param cont_transform A monotonic function to be applied to continuous columns
#' @param cutoff_by When ordinalizing, select cut points by absolute values if `dist` and by quantiles if `quantile`
#' @param qmin Cutoff points are slected in the range of `qmin` and `qmax` quantiles of the data `x`
#' @param qmax See `qmin`.
#' @param num_ord Number of ordinal levels to use
#' @return an ordinal vector.
#' @export
generate_mixed_from_gc <- function(sigma=NULL, n=2000, seed=NULL, var_types = NULL,
                                   cont_transform=NULL, cutoff_by='quantile', qmin=0.05, qmax=0.95, num_ord=5){
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

  X
}

