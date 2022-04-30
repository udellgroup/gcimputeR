#' Truncated normal mean and variance
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean of the normal distribution
#' @param std Standard deviation of the normal distribution
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
moments_truncnorm <- function(mu, std, a, b, tol=1e-6, mean_only=FALSE){
  alpha = (a-mu)/std
  beta = (b-mu)/std
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
    mean_ = mu - r1 * std
    out = list('mean' = mean_)
    #
    if (!mean_only) {
      if (beta >= Inf) r2 = (-alpha * pdf_alpha) / Z
      else if (alpha<=-Inf) r2 = (beta * pdf_beta) / Z
      else r2 = (beta * pdf_beta - alpha * pdf_alpha) / Z
      out[['std']] = std * sqrt(1 - r2 - (r1**2))
    }
  }
  out
}

#' Truncated normal mean and variance
#' @description Compute the mean and variance of one dimensional truncated normal random variable
#' @param mu Mean vector of the normal distribution
#' @param std Standard deviation vector of the normal distribution
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
moments_truncnorm_vec <- function(mu, std, a, b, tol=1e-6, mean_only=FALSE){
  alpha = (a-mu)/std
  beta = (b-mu)/std
  Z = pnorm(beta) - pnorm(alpha)
  p = length(Z)

  if (any(is.infinite(Z)) | min(Z)<0)stop('Invalid input')
  work_loc = which((Z>tol) & (Z<1))
  trivial_loc = Z==1
  fail_loc = Z<=tol

  pdf_beta = dnorm(beta[work_loc])
  pdf_alpha = dnorm(alpha[work_loc])
  if (any(is.infinite(pdf_alpha)) | any(is.infinite(pdf_beta)))stop('Invalid input')
  mean_ = numeric(p)
  r1 = (pdf_beta - pdf_alpha)/Z[work_loc]
  mean_[work_loc] = mu[work_loc] - r1 * std[work_loc]
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

    std_ = numeric(p)
    for (name in names(loc_list)){
      loc = loc_list[[name]]
      abs_loc = work_loc[loc]
      std_[abs_loc] = std[abs_loc] * sqrt(1 - r2_list[[name]] - (r1[loc]**2))
    }
    std_[fail_loc] = Inf
    std_[trivial_loc] = std[trivial_loc]
    out[['std']] = std_
  }

  out
}

#' Compute multivariate truncated normal mean and cov
#'
#' @description  compute multivariate truncated normal mean and cov by sampling
#' @param mean Normal mean vector
#' @param cov Normal covariance matrix
#' @param lower Lower boundary of truncated intervals
#' @param upper Upper boundary of truncated intervals
#' @param n_sample Number of samples to use for sampling methods
#' @return A list containing
#' \describe{
#'   \item{\code{mean}}{Mean vector}
#'   \item{\code{cov}}{Covariance matrix}
#' }
#' @export
get_trunc_2dmoments <- function(mean, cov, lower, upper, method = 'Sampling', n_sample=5000){
  p = length(mean)
  if (p==1){
    out = moments_truncnorm(c(mean), c(cov), c(lower), c(upper))
    return(list(mean = out$mean, cov = matrix(out$std^2, 1, 1)))
  }
  switch (method,
          'Sampling' = {
            z = TruncatedNormal::rtmvnorm(n = n_sample, mu = mean, sigma = cov, lb = lower, ub = upper)
            if (length(dim(z)) != 2 | ncol(z)!=p) stop('unexpected sample dimension')
            r = list('mean' = colMeans(z), 'cov' = cov(z))
          },
          'Diagonal' = {
            r = moments_truncnorm_vec(mu = mean, std = sqrt(diag(cov)), a = lower, b = upper)
            r = list('mean'=r$mean, 'cov'=diag(r$std^2))
          },
          stop(paste0("Invalid method vlaue: ", method))
  )
  return(r)
}

"Explicit = {
            r = tmvtnorm::mtmvnorm(mean = mean, sigma = cov, lower = lower, upper = upper)
            names(r) <- c('mean', 'cov')
          },"
