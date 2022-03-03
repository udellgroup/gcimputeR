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
    #r_lower = qnorm(r_lower * n/(n+1))
    #r_upper = qnorm(r_upper * n/(n+1))
    r_lower = qnorm(r_lower)
    r_upper = qnorm(r_upper)
    return(list(r_lower=r_lower, r_upper=r_upper))
  }
}

#' Transform imputed Z matrix to original scale X
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
  k = length(d_index)
  Xnew = X

  # ordinal
  if (k > 0){
    for (j in 1:k){
      miss_ind = which(is.na(X[,j]))
      # decimal to integer 1 to observed length
      xmis_loc = ceiling(pnorm(Z[miss_ind,j]) * (n-length(miss_ind)))
      # If extremely small value appear in Z[miss_ind, ] such that xmis_loc contains 0,
      # replace it with 1 (impute using the smallest observation)
      xmis_loc[xmis_loc == 0] = 1
      Xnew[miss_ind,j] = sort(X[-miss_ind,j])[xmis_loc]
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



