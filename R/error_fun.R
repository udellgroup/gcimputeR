#' compute the mean absolute error
#'
#' @description compute the mean absolute error on entries which are missing in \code{xobs} but observed in \code{xtrue}
#' @param xhat imputed data matrix or vector
#' @param xobs incomplete observed data matrix or vector
#' @param xtrue complete true data matrix or vector
#' @param round whether round the values of \code{xhat} to integers
#' @return computed error
#' @export
cal_mae = function(xhat, xobs=NULL, xtrue, round = FALSE){
  xhat = as.numeric(as.matrix(xhat))
  xobs = as.numeric(as.matrix(xobs))
  xtrue = as.numeric(as.matrix(xtrue))
  if(round) xhat = round(xhat)
  if (is.null(xobs)) loc = 1:length(xhat) else loc = which(is.na(xobs) & (!is.na(xtrue)))
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
  loc = which(is.na(xobs))
  if (relative) scale = sqrt(mean((xtrue[loc])^2)) else scale = 1
  if (is.null(xobs)) loc = 1:length(xhat) else loc = which(is.na(xobs) & (!is.na(xtrue)))
  sqrt(mean((xhat[loc] - xtrue[loc])^2)) / scale
}

#' compute the scaled mean squared error(SMAE) scaled by the MAE of column
#'
#' @description For each column, compute the MAE of imputed value divided by the MAE of column median imputation
#' @param xhat complete imputed data matrix
#' @param xobs incomplete observed data matrix
#' @param xtrue complete true data matrix
#' @param round whether round the values of \code{xhat} to integers
#' @return a vector with SMAE for each column
#' @export
cal_mae_scaled = function(xhat, xobs, xtrue, round = FALSE){
  n = dim(xtrue)[1]
  p = dim(xtrue)[2]
  xobs = as.numeric(as.matrix(xobs))
  dim(xobs) = c(n,p)
  xhat = as.numeric(as.matrix(xhat))
  dim(xhat) = c(n,p)
  med = apply(xobs,2,median, na.rm=TRUE)
  err.imp = numeric(p)

  for (j in 1:p){
    err.med = cal_mae(xhat = rep(med[j],n), xobs = xobs[,j], xtrue = xtrue[,j])
    err.imp[j] = cal_mae(xhat = xhat[,j], xobs = xobs[,j], xtrue = xtrue[,j], round = round)/err.med
  }
  err.imp
}

