simple_imp <- function(Z, col_mean=FALSE, std=0){
  Zimp = Z
  for (j in 1:ncol(Z)){
    miss = is.na(Z[,j])
    if (col_mean) mu = mean(Z[!miss,j]) else mu =0
    if (std == 0) Zimp[miss,j] = mu
    else Zimp[miss,j] = rnorm(sum(miss), mu, std)
  }
  Zimp
}

#' EM algorithm to fit Gaussian copula
#'
#' @description  fit the Gaussian copula model from incomplete mixed data
#' @param Z Transformed latent matrix
#' @inheritParams initZ_interval_truncated
#' @inheritParams observed_to_latent
#' @inheritParams impute_mixedgc
#' @param corr_min_eigen  If the minimal eigenvalue of a correlation estimate is below \code{corr_min_eigen}, it will be regularized to have minimal eigenvalue equal to \code{corr_min_eigen}
#' @param dcat_index Boolean vector with \code{TRUE} at categorical dimensions
#' @param cat_input Input for categorical dimensions
#' @param start Initial value of copula correlation
#' @param scale_to_corr Whether to scale a covariance into a correlation matrix in each EM iteration. For development purpose. Use with caution.
#' @return A list containing fitted copula correlation matrix, the likelihood(objective function), Z matrix with updated ordinal entries and a complete imputed Z matrix.
#' \describe{
#'   \item{\code{corr}}{Fitted copula correlation matrix}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Z}}{Incomplete \code{Z} with approximated observed ordinal mean}
#'   \item{\code{Zimp}}{Complete \code{Z} with observed entries the same as \code{Zobs} and missing entries imputed}
#' }
#' @export
em_mixedgc = function(Z, Lower, Upper,
                      d_index, dcat_index=NULL,
                      cat_input = NULL,
                      start=NULL,
                      trunc_method='Iterative', n_sample=5000, n_update=1,
                      maxit=50, eps=0.01, verbose=FALSE, runiter=0,
                      corr_min_eigen = 0.01,
                      scale_to_corr=TRUE){
  # Initialize with mean imputation for Z_continuous
  Z_lower = Lower
  Z_upper = Upper
  p = length(d_index)
  c_index = !d_index

  if(is.null(start)){
    Z_meanimp = Z
    Z_meanimp[is.na(Z_meanimp)] = 0
    R = cor(Z_meanimp)
    rm(Z_meanimp)
  }else{
    R = start$R
  }
  o_eigen = eigen(R)
  if (min(o_eigen$values) <0 ){
    message('Bad initialization: Tiny negative eigenvalue potentially due to colinearity')
    values = o_eigen$values
    values[values<corr_min_eigen] = corr_min_eigen
    R = cov2cor(o_eigen$vectors %*% diag(values) %*% t(o_eigen$vectors))
  }
  R = regularize_corr(R, corr_min_eigen = corr_min_eigen, verbose = verbose,
                      prefix = 'Bad initialization potentially due to colinearity: ')
  stopifnot(min(eigen(R)$values)>0)
  if (!is.null(cat_input)) R = project_to_nominal_corr(R, cat_input$cat_index_list)

  Zimp = Z
  l=0
  loglik = NULL

  repeat{
    l = l+1
    est_iter = latent_operation('em',
                                Z, Z_lower, Z_upper,
                                d_index = d_index, dcat_index = dcat_index,
                                cat_input = cat_input,
                                corr = R,
                                trunc_method = trunc_method, n_sample=n_sample, n_update=n_update)
    Z = est_iter$Z
    Zimp = est_iter$Zimp
    R1 = est_iter$corr
    if (scale_to_corr){
      R1 = cov2cor(R1)
      R1 = regularize_corr(R1, corr_min_eigen = corr_min_eigen, verbose = verbose)
      if (!is.null(cat_input)) R1 = project_to_nominal_corr(R1, cat_input$cat_index_list)
      #if (!is.null(cat_input)) R1 = project_to_nominal_corr_simple(R1, cat_input$cat_index_list)
    }
    err = norm(R1-R, type = 'F')/norm(R, type = 'F')
    loglik = c(loglik, est_iter$loglik)
    # update
    R = R1
    if (verbose){
      print(paste0('Iteration ', l, ': ', 'copula parameter change ', round(err, 4), ', likelihood ', round(est_iter$loglik, 4)))
    }
    # determine convergence
    if (runiter==0){
      if (err<eps) break
      if (l > maxit){
        warning('Max iter reached in EM')
        break
      }
    }else{
      if (l>=runiter) break
    }

  }
  return(list(corr=R, loglik=loglik, Z = Z, Zimp = est_iter$Zimp))
}


#' EM algorithm to fit low rank Gaussian copula
#'
#' @description  fit the Gaussian copula model from incomplete mixed data
#' @param rank the number of latent factors
#' @param Z_continuous Continuous columns
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param start Initial value of copula correlation matrix. Default is \code{NULL}.
#' @param maxit maximum number of iterations
#' @param eps Convergence threshold
#' @param verbose If true, output each iteration's detail
#' @return A list containing fitted copula parameters, the likelihood (objective function), Z matrix with updated ordinal entries and the conditional variance corresponding to the observed Z matrix.
#' \describe{
#'   \item{\code{W}}{Fitted latent low rank subspace matrix}
#'   \item{\code{sigma}}{Fitted noise variance}
#'   \item{\code{loglik}}{The log-likelihood achieved during iteration.}
#'   \item{\code{Zobs}}{Incomplete \code{Z} with approximated observed ordinal entries}
#'   \item{\code{C}}{The conditional variance corresponding to the observed Z matrix. Useful for quantifying imputation uncertainty.}
#'   \item{\code{S}}{Required quantity to impute the Z matrix}
#' }
#' @author Yuxuan Zhao, \email{yz2295@cornell.edu} and Madeleine Udell, \email{udell@cornell.edu}
#' @references Zhao, Y., & Udell, M. (2020). Matrix Completion with Quantified Uncertainty through Low Rank Gaussian Copula. arXiv preprint arXiv:2006.10829.
#' @export
em_mixedgc_ppca = function(rank, Z_continuous, r_lower, r_upper,
                           start =NULL, maxit=100, eps=0.01, verbose = FALSE){
  if (is.null(Z_continuous)){
    p = dim(r_upper)[2]
    k = p
  }
  else{
    if (is.null(r_upper)){
      p = dim(Z_continuous)[2]
      k=0
    }
    else{
      p = dim(r_upper)[2] +  dim(Z_continuous)[2]
      k = dim(r_upper)[2]
    }
  }

  # Initialize observed ordinal values
  if (k > 0){
    Z_ordinal = initZ_interval_truncated(Lower = r_lower, Upper = r_upper)
  } else Z_ordinal = NULL

  # Initialize with SVD imputation for Z
  if(is.null(start)){
    Z_meanimp = cbind(Z_ordinal, Z_continuous)
    Z_meanimp[is.na(Z_meanimp)] = 0

    # approximate Z_meanimp by a low rank SVD
    # and correspondingly update Z_ordinal
    Z_meanimp = impute_init(Z_meanimp, rank, r_upper, r_lower)
    ind = which(!is.na(Z_ordinal))
    Z_ordinal[ind] = Z_meanimp[ind]
    Z = cbind(Z_ordinal, Z_continuous)
    Z_meanimp[!is.na(Z)] = Z[!is.na(Z)]

    R = cor(Z_meanimp)
    est = svd(R)
    rm(Z_meanimp, R)

    sigma = mean(est$d[-(1:rank)])
    W = est$u[,1:rank] %*% diag(sqrt(est$d[1:rank] - sigma))

    # scaling
    est = scale_corr(W,sigma)
    W = est$W
    sigma = est$sigma
    rm(est)
  }else{
    W = start$W
    sigma = start$sigma
  }

  Z = cbind(Z_ordinal, Z_continuous)
  l=0
  loglik = NULL
  #diffW = list()
  #estvar = NULL
  repeat{
    l = l+1
    est_iter = em_mixedgc_ppca_iter(Z, r_lower, r_upper, W, sigma)
    Z = est_iter$Zobs
    W1 = est_iter$W
    loglik = c(loglik, est_iter$loglik)
    #
    err = sqrt(sum((W1-W)^2)/sum(W^2))
    #diffW[[l]] = c(err, grassman_dist(W1,W))
    #estvar = c(estvar, est_iter$sigma)
    if (verbose){
      print(paste('update error:', err, ' estimated var:',est_iter$sigma))
      print(paste('likelihood: ', est_iter$loglik))
    }
    if (err<eps) break
    if (l > maxit){
      warning('Max iter reached in EM')
      break
    }
    W = W1
    sigma = est_iter$sigma
  }

  return(list(W=W, sigma = sigma, loglik=loglik, Zobs = Z, S = est_iter$S, C = est_iter$C))
}
