#' Generate values for latent \code{Z} at observed ordinal entires
#'
#' @description  Each observed ordinal entry of \code{Z} follows truncated normal distribution with mean 0, variance 1 and truncated interval provided in \code{lower} and \code{upper}. Fill out those entries by their univariate mean if \code{method = "univariate_mean"} and by random sampples if \code{method = "sampling"}
#' @param Lower Lower boundary of truncated intervals
#' @param Upper Upper boundary of truncated intervals
#' @param seed random seed used
#' @param method method for initializing the mean
#' @return Matrix \code{Z} with the same shape as \code{lower}
#' @export
initZ_ordinal = function(Lower, Upper, seed = NULL, method = 'univariate_mean'){
  if (!is.null(seed)) set.seed(seed)
  # input: truncated interval, mean and variance
  # output: for each entry, a random number following the truncated normal with mean 0 and var 1 is returned
  # complexity: the number of observed ordinal entries
  r_upper = Upper
  r_lower = Lower
  n = dim(r_upper)[1]
  k = dim(r_upper)[2]
  Z = matrix(NA, n, k)

  obs_indices = !is.na(r_lower)
  u_lower = pnorm(r_lower[obs_indices])
  u_upper = pnorm(r_upper[obs_indices])

  if (min(u_upper-u_lower)<=0){
    loc = which.min(u_upper-u_lower)
    print(paste('Min of upper - lower', u_upper[loc]-u_lower[loc]))
    print(paste('where upper is', u_upper[loc], 'and lower is', u_lower[loc]))
    stop()
  }
  if (min(u_lower)<0)stop(paste('Invalid min of lower', min(u_lower)))
  if (max(u_upper)>1)stop(paste('Invalid max of upper', max(u_upper)))

  switch (method,
          'sampling' = {
            Z[obs_indices] = qnorm(purrr::map2_dbl(u_lower, u_upper, runif, n=1))
          },
          'univariate_mean' = {
            l = sum(obs_indices)
            out = moments_truncnorm_vec(mu = numeric(l),std = 1+numeric(l),
                                        a = r_lower[obs_indices], b = r_upper[obs_indices], mean_only = TRUE)
            Z[obs_indices] = out$mean
          },
          stop('invalid method')
  )
  Z
}
