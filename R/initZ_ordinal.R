#' Simulate random value for observed ordinal entires of \code{Z}
#'
#' @description  Each observed ordinal entry of \code{Z} follows truncated normal distribution with mean 0, variance 1 and truncated interval provided in \code{r_lower} and \code{r_upper}. Exact values of observed ordinal entries of \code{Z} are needed to implement the EM algorith. The initial values are ramdom values from truncated normal distribution.
#' @param r_lower Lower boundary of truncated intervals for ordinal columns
#' @param r_upper Upper boundary of truncated intervals for ordinal columns
#' @param seed random seed used. Default is 1
#' @return A matrix containing ordinal columns of \code{Z} whose entries are randomly generated but lie in the truncated interval.
#' @export
initZ_ordinal = function(r_lower, r_upper, seed = 1){
  # input: truncated interval, mean and variance
  # output: for each entry, a random number following the truncated normal with mean 0 and var 1 is returned
  # complexity: the number of observed ordinal entries
  n = dim(r_upper)[1]
  k = dim(r_upper)[2]
  Z = r_upper
  set.seed(seed)
  for (i in 1:n){
    i_obs = which(!is.na(Z[i,]))

    # only implement when observing ordinal values
    if (length(i_obs) > 0){
      u_lower = pnorm(r_lower[i,i_obs])
      u_upper = pnorm(r_upper[i,i_obs])
      for (j in 1:length(i_obs)){
        if ((u_lower[j] < 1) & (u_upper[j] > 0)){
          u = runif(1, u_lower[j], u_upper[j])
          if (is.na(u)) stop(paste('NAs produced in simulating u: '), u_lower[j], u_upper[j])
          Z[i,i_obs[j]] = qnorm(u)
        }
      }
    }
  }
  Z
}
