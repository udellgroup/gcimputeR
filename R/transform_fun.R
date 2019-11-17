#' Transform original observation X to copula observation Z
#'
#' @description  For each continuous column of \code{X}, compute its corresonding value of \code{Z}; For each ordinal column of \code{Z}, compute its corresponding truncated interval of \code{Z}.
#' @param X original observation matrix
#' @param type Speficy data type. Either \code{continuous} or \code{ordinal}.
#' @return A list containing
#' \describe{
#'   \item{\code{r_val}}{Return when \code{type=continuous}. Corresponding value}
#'   \item{\code{r_lower}}{Return when \code{type=ordinal}. Lower boundary of corresponding truncated interval}
#'   \item{\code{r_upper}}{Return when \code{type=ordinal}. Upper boundary of corresponding truncated interval}
#' }
#' @export
range_transform = function(X, type){
  n = dim(X)[1]
  p = dim(X)[2]

  # for continuous dimensions, return a transformed value for each observation
  if (type == 'continuous'){
    r_val = X
    for (j in 1:p){
      fun = ecdf(X[,j])
      r_val[,j] = fun(X[,j])
    }
    r_val = qnorm(r_val * n/(n+1))
    return(list(r_val = r_val))
  }

  # for ordinal dimensions, return a transformed interval for each observation
  if (type == 'ordinal'){
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
    r_lower = qnorm(r_lower * n/(n+1))
    r_upper = qnorm(r_upper * n/(n+1))
    return(list(r_lower=r_lower, r_upper=r_upper))
  }
}

#' Transform imputed Z matrix to original scale X
#'
#' @description  Impute orginal incomplete observation using imputed \code{Z} matrix and marginal distribution information drawn from original observation \code{X}
#' @param Z imputed Z matrix
#' @param X original observation matrix
#' @param r_upper Speficy data type. Either \code{continuous} or \code{ordinal}.
#' @param d_index ordinal dimension indexes.
#' @return Imputed data matrix
#' @export
Ximp_transform = function(Z, X, r_upper, d_index){
  # for continuous, use default empirical quantile
  # deal with smaller dimensional r_upper
  n = dim(Z)[1]
  p = dim(Z)[2]
  k = length(d_index)
  Xnew = X

  # ordinal
  if (k > 0){
    for (j in 1:k){
      miss_ind = which(is.na(X[,j]))
      for (i in miss_ind){
        indset = which(r_upper[,j] >= Z[i,j])
        if (length(indset) == 0) Xnew[i,j] = max(X[,j], na.rm = TRUE)
        else{
          ind = indset[which.min(r_upper[indset,j])]
          Xnew[i,j] = X[ind,j]
        }
      }
    }
  }

  # continuous
  if (k < p){
    for (j in (k+1):p){
      miss_ind = which(is.na(X[,j]))
      Xnew[miss_ind,j] = quantile(X[,j], pnorm(Z[miss_ind,j]), na.rm = TRUE)
    }
  }

  Xnew
}

#' Transform a continuous vector to ordinal vector through cutoff function
#'
#' @description  Discretize continuous x to ordinal z with \code{k} levels. The cutoffs are either given or randomly generated.
#' @param x continuous vector
#' @param k number of ordinal levels. Generated ordinal vector takes value from \code{0,\ldots,k-1}
#' @param cutoff cutoffs. Setting it as \code{NULL} will generate random cutoffs.
#' @return an ordinal vector.
#' @export
continuous2ordinal = function(x, k = 2, cutoff = NULL){
  # if cutoff is not provided, a random cutoff sequence will be generated.
  if (k == 2){
    if (is.null(cutoff)){
      repeat{
        cutoff = sample(x,1)
        q = quantile(x, c(0.1,0.9))
        if ((cutoff>q[1]) & (cutoff<q[2])) break # make sure each class has at least 10% data
      }
    }
    x = (x >= cutoff)
  }
  else{
    if (is.null(cutoff)){
      cutoff = seq(min(x)-0.1 *sd(x), max(x) + 0.1*sd(x), length.out = k+1)
    }
    x = cut(x, cutoff, labels = FALSE, include.lowest = TRUE)
  }
  x
}
