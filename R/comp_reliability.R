#' Compute reliability
#'
#' @description Compute reliability for continuous matrices
#' @param X Incomplete data
#' @param fit Fitted LRGC model on data \code{X}
#' @param alpha The confidence leveled used in computing reliability
#' @param type Default \code{1} corresponds to the one used in our paper
#' @export
reliability_cont = function(X, fit, alpha = 0.95, type = 1){
  ct = ct_impute(X, fit, alpha = alpha)

  loc = which(!is.na(ct$upper))
  d =  ct$upper[loc] - ct$lower[loc]
  ximp = fit$Ximp[loc]

  if (type == 1){
    r = numeric(length(loc))
    sumx = sum(ximp^2)
    sumd = sum(d^2)
    for (i in 1:length(loc)) r[i] = (sumd - d[i]^2)/(sumx - ximp[i]^2)
  }
  if (type == 2){
    r = d^2/ximp^2
  }
  if (type == 3){
    r = d^2
  }

  R = array(dim = dim(X))
  R[loc] = r
  R
}


#' Compute reliability
#'
#' @description Compute reliability for ordinal matrices
#' @param X Incomplete data
#' @param fit Fitted LRGC model on data \code{X}
#' @param cuts.fix If \code{NULL}, use the cut points estimated from \code{fit}
#' @export
reliability_ord = function(X, fit, cuts.fix = NULL){
  # output: probability for predicted ordinal entries
  # fit is a returned list from function impute_mixedgc_ppca
  # used objects: W, sigma, cutoffs
  # C (the conditional variance of observed ordinal), Zimp (imputed values of Z),
  # The original observation X is only used to detect the missing locations
  n = dim(X)[1]
  p = dim(X)[2]
  obj = svd(fit$W)
  U = obj$u
  d = obj$d
  rank = length(d)

  sig = fit$sigma
  C = fit$C
  Zimp = fit$Zimp
  cutoffs = fit$cutoffs

  prob = array(NA, dim = c(n,p))

  # compute conditional variance
  for (i in 1:n){
    index_m = which(is.na(X[i,]))
    index_o = which(!is.na(X[i,]))
    Uiobs = matrix(U[index_o,], ncol = rank)
    Uimis = matrix(U[index_m,], ncol = rank)

    #uCu = t(Uiobs) %*% diag(C[i,index_o],nrow = length(index_o)) %*% Uiobs
    uCu = t(Uiobs) %*% (C[i,index_o] * Uiobs)
    dUmis = solve(sig * diag(d^{-2}) + t(Uiobs) %*% Uiobs, t(Uimis))
    # compute variance i.e. diagonal elements of conditional covariance matrix
    for (l in 1:length(index_m)){
      j = index_m[l]
      du = dUmis[,l]
      C[i,j] = sig + sig * sum(du * U[j,]) + sum(du * (uCu %*% du))
    }
  }

  # compute probability
  for (j in 1:p){
    cuts = cutoffs[[j]]
    index_m = which(is.na(X[,j]))
    for (i in index_m){
      if (is.null(cuts.fix)){
        t = min(abs(Zimp[i,j]- cuts))
      }else{
        t = abs(Zimp[i,j]- sort(cuts,decreasing = TRUE)[2])
      }
      prob[i,j] = 1 - C[i,j]/t^2
    }
  }

  prob
}
