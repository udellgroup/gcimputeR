#' Truncated normal mean and variance
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

#' @rdname mean_tnorm
#' @inheritParams mean_tnorm
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

