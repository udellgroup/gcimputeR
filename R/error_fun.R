
#' compute the mean absolute error
#'
#' @description compute the mean absolute error on entries which are missing in \code{xobs} but observed in \code{xtrue}
#' @param xhat imputed data matrix or vector
#' @param xobs incomplete observed data matrix or vector
#' @param xtrue complete true data matrix or vector
#' @param round whether round the values of \code{xhat} to observed integers
#' @return computed error
#' @export
cal_mae = function(xhat, xobs=NULL, xtrue, round = FALSE){
  xhat = as.numeric(as.matrix(xhat))
  xobs = as.numeric(as.matrix(xobs))
  xtrue = as.numeric(as.matrix(xtrue))
  if(round) xhat = round(xhat)
  if (is.null(xobs)) loc = !is.na(xtrue) else loc = which(is.na(xobs) & (!is.na(xtrue)))
  mean(abs(xhat[loc] - xtrue[loc]))
}

#' compute the root mean squared error
#'
#' @description compute the root mean squared error on entries which are missing in \code{xobs} but observed in \code{xtrue}
#' @param xhat imputed data matrix or vector
#' @param xobs incomplete observed data matrix or vector
#' @param xtrue complete true data matrix or vector
#' @param relative whether divide the computed error by the Frobenius norm of \code{xtrue}
#' @return computed error
#' @export
cal_rmse = function(xhat, xobs=NULL, xtrue, relative = TRUE){
  xhat = as.numeric(as.matrix(xhat))
  xobs = as.numeric(as.matrix(xobs))
  xtrue = as.numeric(as.matrix(xtrue))
  if (is.null(xobs)) loc = !is.na(xtrue) else loc = which(is.na(xobs) & (!is.na(xtrue)))
  if (relative) scale = sqrt(mean((xtrue[loc])^2)) else scale = 1
  sqrt(mean((xhat[loc] - xtrue[loc])^2)) / scale
}

#' compute the scaled mean squared error(SMAE) scaled by the MAE of column
#'
#' @description For each column, compute the MAE of imputed value divided by the MAE of column median imputation
#' @param xhat complete imputed data matrix
#' @param xobs incomplete observed data matrix
#' @param xtrue complete true data matrix
#' @param round whether round the values of \code{xhat} to integers
#' @param base_from_true If \code{TRUE}, form baseline imputation from \code{xtrue}. Else, form baseline imputation from \code{xobs}.
#' @param verbose If \code{TRUE}, throw out a warning when perfect baseline imputation appears
#' @return a vector with SMAE for each column
#' @export
cal_mae_scaled = function(xhat, xobs, xtrue,
                          round = FALSE, reduce = TRUE, base_from_true = FALSE, verbose = TRUE){
  scale_off = FALSE
  n = dim(xtrue)[1]
  p = dim(xtrue)[2]
  xobs = to_numeric_matrix(xobs)
  xhat = to_numeric_matrix(xhat)
  xtrue = to_numeric_matrix(xtrue)
  if (base_from_true) xbase = xtrue else xbase = xobs
  med = apply(xbase,2,median, na.rm=TRUE)
  err = numeric(p)

  for (j in 1:p){
    if (!scale_off){
      if (round) medimp = to_nearest_ord(med[j], xobs[,j]) else medimp = med[j]
      err_med = cal_mae(xhat = rep(medimp,n), xobs = xobs[,j], xtrue = xtrue[,j])
    }else err_med = 1
    if (round) ximp = to_nearest_ord_vec(xhat[,j], xobs[,j]) else ximp = xhat[,j]
    err[j] = cal_mae(xhat = ximp, xobs = xobs[,j], xtrue = xtrue[,j])/err_med
  }
  keep = is.finite(err)
  remove = !keep
  if (any(remove) & verbose){
    #stop('Perfect median imputation appears')
    warning(paste0('Perfect majority vote imputation appears in ', sum(remove), ' vars'))
  }
  if (reduce) err = mean(err[keep])
  err
}

#' compute mis-classification rate for categorical data
#'
#' @description compute mis-classification rate for categorical data
#' @param xhat complete imputed data matrix
#' @param xobs incomplete observed data matrix
#' @param xtrue complete true data matrix
#' @export
cal_misclass = function(xhat, xobs, xtrue){
  xobs = as.numeric(as.matrix(xobs))
  xhat = as.numeric(as.matrix(xhat))
  xtrue = as.numeric(as.matrix(xtrue))
  if (is.null(xobs)) loc = !is.na(xtrue) else loc = is.na(xobs) & (!is.na(xtrue))
  1 - mean(xhat[loc] == xtrue[loc])
}

cal_misclass_scaled = function(xhat, xobs, xtrue, base_from_true = FALSE, reduce =TRUE){
  xobs = to_numeric_matrix(xobs)
  xhat = to_numeric_matrix(xhat)
  xtrue = to_numeric_matrix(xtrue)

  if (base_from_true) xbase = xtrue else xbase = xobs
  base = apply(xbase, 2, function(x) names(which.max(table(x))))
  base = as.integer(base)
  p = ncol(xhat)
  err = numeric(p)
  n = nrow(xobs)

  for (j in 1:p){
    err_base = cal_misclass(xhat = rep(base[j],n), xobs = xobs[,j], xtrue = xtrue[,j])
    err[j] = cal_misclass(xhat = xhat[,j], xobs = xobs[,j], xtrue = xtrue[,j])/err_base
  }

  remove = is.infinite(err)
  if (any(remove)){
    #stop('Perfect median imputation appears')
    warning(paste0('Perfect majority vote imputation appears in ', sum(remove), ' vars'))
  }
  if (reduce) err = mean(err[!remove])
  err
}
