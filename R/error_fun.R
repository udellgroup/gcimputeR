#' compute mean absolute error
#'
#' @param xhat imputed value vector
#' @param xtrue true value vector
#' @param round whether round the values of \code{xhat} to integers
#' @return computed error
#' @export
cal_mae = function(xhat,xtrue, round = FALSE){
  xhat = as.numeric(xhat)
  if(round) xhat = round(xhat)
  mean(abs(xhat - xtrue))
}

#' compute root mean squared error
#'
#' @param xhat imputed value vector
#' @param xtrue true value vector
#' @return computed error
#' @export
cal_rmse = function(xhat,xtrue){
  xhat = as.numeric(xhat)
  sqrt(mean((xhat - xtrue)^2))
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
    ind = which(is.na(xobs[,j]) & (!is.na(xtrue[,j]))) # test points in dimension j
    err.med = cal_mae(med[j], xtrue[ind,j])
    err.imp[j] = cal_mae(xhat[ind,j], xtrue[ind,j], round = round)/err.med
  }
  err.imp
}

