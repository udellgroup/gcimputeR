#' Compute imputation error
#'
#' @description Compute the imputation error on entries which are missing in \code{xobs} but observed in \code{xtrue}
#' @param xhat imputed data matrix or vector
#' @param xobs incomplete observed data matrix or vector
#' @param xtrue complete true data matrix or vector
#' @param round whether to round the values of \code{xhat} to observed integers in \code{xobs}
#' @param relative Only used for \code{cal_rmse}. If \code{TRUE}, normalize the computed error by the Frobenius norm of \code{xtrue}
#' @param reduce Only used for \code{cal_smae}. Return the average error if \code{TRUE} else the whole error vector
#' @param verbose Only used for \code{cal_smae}. If \code{TRUE}, throw out a warning when perfect baseline imputation appears
#' @return Computed imputation error
#' @name comp_error
NULL
#> NULL

#' @describeIn comp_error Compute the mean absolute error
#' @export
cal_mae = function(xhat, xobs, xtrue, round = FALSE){
  xhat = as.numeric(as.matrix(xhat))
  xobs = as.numeric(as.matrix(xobs))
  xtrue = as.numeric(as.matrix(xtrue))
  if(round) xhat = round(xhat)
  #if (is.null(xobs)) loc = !is.na(xtrue) else loc = which(is.na(xobs) & (!is.na(xtrue)))
  loc = which(is.na(xobs) & (!is.na(xtrue)))
  mean(abs(xhat[loc] - xtrue[loc]))
}

#' @describeIn comp_error Compute the root mean squared error
#' @export
cal_rmse = function(xhat, xobs, xtrue, relative = TRUE){
  xhat = as.numeric(as.matrix(xhat))
  xobs = as.numeric(as.matrix(xobs))
  xtrue = as.numeric(as.matrix(xtrue))
  #if (is.null(xobs)) loc = !is.na(xtrue) else loc = which(is.na(xobs) & (!is.na(xtrue)))
  loc = which(is.na(xobs) & (!is.na(xtrue)))
  if (relative) scale = sqrt(mean((xtrue[loc])^2)) else scale = 1
  sqrt(mean((xhat[loc] - xtrue[loc])^2)) / scale
}

#' Wrapper for cal_smae
#' @export
#' @keywords internal
cal_mae_scaled = function(...) cal_smae(...)

#' @describeIn comp_error Compute the MAE scaled by median imputation MAE (at each column)
#' @export
cal_smae = function(xhat, xobs, xtrue,
                          round = FALSE, reduce = TRUE, verbose = TRUE){
  # If \code{TRUE}, form baseline imputation from \code{xtrue}. Else, form baseline imputation from \code{xobs}.
  base_from_true = FALSE
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

#' @describeIn comp_error Compute the mis-classification rate for categorical data
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
