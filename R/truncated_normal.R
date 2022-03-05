#' Truncated normal mean
#'
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param a Left endpoint of the truncated interval
#' @param b Right enpoint of the truncated interval
#' @export
mean_tnorm = function(mu,sigma,a,b){
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  Z = pnorm(beta) - pnorm(alpha)
  mu + (dnorm(alpha) - dnorm(beta)) / Z * sigma
}

#' Truncated normal variance
#'
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param a Left endpoint of the truncated interval
#' @param b Right enpoint of the truncated interval
#' @export
var_tnorm = function(mu,sigma,a,b){
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  Z = pnorm(beta) - pnorm(alpha)
  if(a == -Inf) return(sigma^2 * (1 - beta * dnorm(beta)/Z - (dnorm(beta)/Z)^2))
  if(b == Inf) return(sigma^2 * (1 + alpha * dnorm(alpha)/Z - (dnorm(alpha)/Z)^2))
  if(is.finite(a) & is.finite(b)){
    return(sigma^2 * (1 + (alpha*dnorm(alpha) - beta*dnorm(beta))/Z - ((dnorm(alpha) - dnorm(beta))/Z)^2))
  }
}

#' Truncated normal mean and variance
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param a Left endpoint of the truncated interval
#' @param b Right endpoint of the truncated interval
#' @param tol Computation is only conducted if the truncated normal probability is larger than \code{tol}
#' @param mean_only Only compute the truncated normal variance when `mean_only=FALSE`
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{The truncated normal mean}
#'   \item{\code{var}}{The truncated normal varaince or `NULL` }
#' }
#' @export
moments_truncnorm <- function(mu, sigma, a, b, tol=1e-6, mean_only=FALSE){
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  Z = pnorm(beta) - pnorm(alpha)

  if (is.infinite(Z))stop('Invalid input')

  if (Z<tol){
    out = list('mean'=Inf, 'var'=Inf)
  }else{
    pdf_beta = dnorm(beta)
    pdf_alpha = dnorm(alpha)
    if (is.infinite(pdf_alpha) | is.infinite(pdf_beta)) stop('Invalid input')
    if (is.infinite(alpha) & is.infinite(beta)) stop('Invalid input')
    r1 = (pdf_beta - pdf_alpha)/Z
    mean_ = mu - r1 * sigma
    out = list('mean' = mean_)
    #
    if (!mean_only) {
      if (beta >= Inf) r2 = (-alpha * pdf_alpha) / Z
      else if (alpha<=-Inf) r2 = (beta * pdf_beta) / Z
      else r2 = (beta * pdf_beta - alpha * pdf_alpha) / Z
      var_ = (sigma**2) * (1 - r2 - (r1**2))
      out[['var']] = var_
    }
  }
  out
}

#' Truncated normal mean and variance
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean vector of the normal distribution
#' @param sigma Standard deviation vector of the normal distribution
#' @param a Left endpoint vector of the truncated interval
#' @param b Right endpoint vector of the truncated interval
#' @param tol Computation is only conducted if the truncated normal probability is larger than \code{tol}
#' @param mean_only Only compute the truncated normal variance when `mean_only=FALSE`
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{The truncated normal mean with the same length of `mu` }
#'   \item{\code{var}}{The truncated normal varaince with the same length of `mu` or `NULL` }
#' }
#' @export
moments_truncnorm_vec <- function(mu, sigma, a, b, tol=1e-6, mean_only=FALSE){
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  Z = pnorm(beta) - pnorm(alpha)
  p = length(Z)

  if (any(is.infinite(Z)) | min(Z)<0)stop('Invalid input')
  work_loc = (Z>tol) & (Z<1)
  trivial_loc = Z==1
  fail_loc = Z<=tol

  pdf_beta = dnorm(beta[work_loc])
  pdf_alpha = dnorm(alpha[work_loc])
  if (any(is.infinite(pdf_alpha)) | any(is.infinite(pdf_beta)))stop('Invalid input')
  mean_ = numeric(p)
  r1 = (pdf_beta - pdf_alpha)/Z[work_loc]
  mean_[work_loc] = mu[work_loc] - r1 * sigma[work_loc]
  mean_[fail_loc] = Inf
  mean_[trivial_loc] = mu[trivial_loc]
  out = list('mean' = mean_)

  if (!mean_only){
    loc_list = list()
    r2_list = list()

    beta_work = beta[work_loc]
    alpha_work = alpha[work_loc]
    Z_work = Z[work_loc]
    #
    loc = beta_work >= Inf
    if (any(loc)){
      r2_list[['inf_beta']] = (-alpha_work[loc] * pdf_alpha[loc]) / Z_work[loc]
      loc_list[['inf_beta']] = loc
    }
    #
    loc = alpha_work <= -Inf
    if (any(loc)){
      r2_list[['inf_alpha']] = (beta_work[loc] * pdf_beta[loc]) / Z_work[loc]
      loc_list[['inf_alpha']] = loc
    }
    #
    loc = (beta_work < Inf) & (alpha_work > -Inf)
    if (any(loc)){
      r2_list[['finite']] = (beta_work[loc] * pdf_beta[loc] - alpha_work[loc] * pdf_alpha[loc]) / Z_work[loc]
      loc_list[['finite']] = loc
    }

    var_ = numeric(p)
    for (name in names(loc_list)){
      loc = loc_list$name
      abs_loc = work_loc[loc]
      var_[abs_loc] = (sigma[abs_loc]**2) * (1 - r2_list$name - (r1[loc]**2))
    }
    var_[fail_loc] = Inf
    var_[trivial_loc] = sigma[trivial_loc]**2
    out[['var']] = var_
  }

  out
}
