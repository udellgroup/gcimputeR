#' Compute categorical frequency
#'
#' @description  For a  categorical matrix, compute the empirical frequency for each variable
#' @param X_cat An integer encoded categorical matrix
#' @export
get_cat_index_freq <- function(X_cat){
  p_cat = ncol(X_cat)
  freq = vector('list', p_cat)
  nlevel = numeric(p_cat)
  for (j in 1:p_cat){
    #count = tabulate(X_cat[,j])
    count = unname(table(X_cat[,j]))
    freq[[j]] = count/sum(count)
    nlevel[j] = length(count)
  }
  list(freq = freq, nlevel = nlevel)
}

softmax_colwise <- function(Z){
  Zexp = exp(Z - apply(Z,1,max))
  Zexp_rsum = rowSums(Zexp)
  r = Zexp/rep(Zexp_rsum,ncol(Z))
  r
}

E_softmax_MC <- function(mu, beta, n_MC=2000, seed=NULL, old=FALSE){
  if (old) return(E_softmax_MC_old(mu, beta, n_MC, seed))
  mu = c(0, mu)
  if (!is.null(seed)) set.seed(seed)
  d = length(mu)
  Z = matrix(rnorm(n_MC * d), n_MC) + rep(mu, each = n_MC)
  pi_x = softmax_colwise(Z*beta)
  val = colMeans(pi_x)
  Inkk = array(0, dim = c(n_MC, d, d))
  pi_x_1st = array(0, dim = c(n_MC, d, d))
  pi_x_2nd = array(0, dim = c(n_MC, d, d))
  for (i in 1:n_MC) Inkk[i,,] = diag(d)
  for (i in 1:d){
    pi_x_1st[,i,] = pi_x
    pi_x_2nd[,,i] = pi_x
  }
  jac = (Inkk - pi_x_1st) * pi_x_2nd
  jac = apply(jac, c(2,3), mean) * beta
  list(val = val, jac = jac)
}

E_softmax_MC_old <- function(mu, beta, n_MC=2000, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  d = length(mu)
  Z = matrix(0, n_MC, d+1)
  Z[,-1] = rnorm(n_MC * d) + rep(mu, each = n_MC)
  #Z = matrix(rnorm(n_MC * (d+1)), n_MC) + rep(c(0,mu), each = n_MC)
  pi_x = softmax_colwise(Z*beta)
  val = colMeans(pi_x)
  Inkk = array(0, dim = c(n_MC, d+1, d+1))
  pi_x_1st = array(0, dim = c(n_MC, d+1, d+1))
  pi_x_2nd = array(0, dim = c(n_MC, d+1, d+1))
  for (i in 1:n_MC) Inkk[i,,] = diag(d+1)
  for (i in 1:(d+1)){
    pi_x_1st[,i,] = pi_x
    pi_x_2nd[,,i] = pi_x
  }
  jac = (Inkk - pi_x_1st) * pi_x_2nd
  jac = apply(jac, c(2,3), mean) * beta
  list(val = val, jac = jac)
}

get_solve_mu <- function(prob, beta, n_MC=2000, seed=11, old=FALSE){
  if (!is.null(seed)) set.seed(seed)
  force(prob)
  force(beta)
  force(n_MC)
  force(seed)
  force(old)
  f_val <- function(mu){
    out = E_softmax_MC(mu, beta, n_MC = n_MC, seed = seed, old = old)
    val = out$val
    val[-1] - prob[-1]
  }
  f_jac <- function(mu){
    out = E_softmax_MC(mu, beta, n_MC = n_MC, seed = seed, old = old)
    jac = out$jac
    if (old) jac = jac[-1,-1]
  }
  list(f_val = f_val, f_jac = f_jac)
}

init_mu <- function(prob, old=FALSE){
  d = length(prob)
  prob = prob[1] / (prob[-1] + prob[1])
  if (old){
    -qnorm(prob)
  }else{
    -sqrt(2)*qnorm(prob)
  }
}

solve_nominal_mu <- function(prob, beta = 1000, n_MC = 10000, seed = 101, inits = NULL, eps = 1e-4, old=FALSE){
  n_MC = c(5000, 10000, 15000, 20000, 25000, 30000)
  best_precis = 100
  best_r = NULL

  mu0 = init_mu(prob, old=old)
  d = length(mu0)

  if (is.null(inits)){
    p_seq = 0:5
    l = length(p_seq)
    inits = matrix(mu0, nrow = l, ncol = d, byrow = TRUE)
    inits = inits / matrix(rep(2^p_seq, each = d), nrow = l, byrow = TRUE)
    inits[l,] = 0
  }
  l = nrow(inits)

  for (m in n_MC){
    # using m MC samples
    flist = get_solve_mu(prob = prob, beta = beta, n_MC = m, seed=seed, old=old)
    # from previous fit
    if (best_precis < 0.1){
      r = rootSolve::multiroot(f= flist$f_val, start = best_r$root, jacfunc = flist$f_jac, verbose=FALSE)
      ifbest = r$estim.precis < best_precis
      stopifnot(!is.na(ifbest))
      if (ifbest){
        best_precis = r$estim.precis
        best_r = r
      }
    }
    if (best_precis < eps) break

    # solving from simple starting points
    for (i in 1:l){
      r = rootSolve::multiroot(f= flist$f_val, start = inits[i,], jacfunc = flist$f_jac, verbose=FALSE)
      ifbest = r$estim.precis < best_precis
      stopifnot(!is.na(ifbest))
      if (ifbest){
        best_precis = r$estim.precis
        best_r = r
      }
      if (best_precis < eps) break
    }
    if (best_precis < eps) break
  }
  best_r
}

#' Categorical marginal estimation
#'
#' @description  Fit the mean of Gaussian-Max distribution for given frequency
#' @param freq_list A list of category frequency
#' @param beta Constant in softmax computation
#' @param n_MC Number of samples used to approximate the expectation of softmax
#' @param seed random seed
#' @param eps The largest accepted marginal estimation error
#' @param verbose Whether to print progress information
#' @param old Use previous formulation?
#' @export
get_cat_mu <- function(freq_list, beta = 1000, n_MC = 5000, seed = 101,  eps = 1e-4,
                       verbose = FALSE, old = FALSE){
  inits = NULL
  quiet_solve=TRUE
  if (verbose) print('Starting categorical marginal estimation: ')
  d = sum_list_len(freq_list)
  lf = length(freq_list)
  if (old) mu = numeric(d-lf) else mu = numeric(d)
  start = 0
  if (old) fs = numeric(d-lf) else fs = numeric(d)
  est_precis = numeric(lf)
  if (quiet_solve) f = purrr::quietly(solve_nominal_mu)
  else f = solve_nominal_mu
  for (i in seq_along(freq_list)){
    out = f(freq_list[[i]],
            beta = beta, n_MC = n_MC, seed = seed, inits = inits, eps = eps)
    if (quiet_solve) ri = out$result else ri = out
    mui = ri$root
    if (old) index = start + (1:length(mui)) else index = start + 1 + (1:length(mui))
    mu[index] = mui
    fs[index] = ri$f.root
    est_precis[i] = ri$estim.precis
    if (old) start = start + length(mui) else start = start + 1 + length(mui)
    if (verbose) print(paste0(i, ' of ', lf, ' categorical marginal estimations finished'))
  }
  if (any(est_precis > eps)) warning('Some nominal mean estimation does not reach desired precision')
  list(mu = mu, f_root = fs, est_precis = est_precis)
}

"
p0 = c(0.5, 0.2 ,0.3)
for (old in c(TRUE, FALSE)){
  flist = get_solve_mu(prob = p0, beta = 1000, n_MC = 10000, seed=101, old=old)
  mu0 = init_mu(p0, old=old)
  r = rootSolve::multiroot(f= flist$f_val, start = mu0, jacfunc = flist$f_jac)
  print(paste0('OLD paradigm: ', old))
  print(r)
}
"

'
[1] "OLD paradigm: TRUE"
$root
[1] -0.6888746 -0.4310904
$f.root
[1] 2.525551e-07 2.384466e-07
$iter
[1] 5
$estim.precis
[1] 2.455009e-07

[1] "OLD paradigm: FALSE"
$root
[1] -0.7434508 -0.4236463
$f.root
[1] -1.179909e-11  2.229139e-10
$iter
[1] 5
$estim.precis
[1] 1.173565e-10
'

"
source('~/gcimputeR/Development_tests/nominal_develop/func_data_generation.R')
d_gen = gen_nominal_copula(n=1000, p_cat_vec = c(3, 4, 5), p_noncat = 0, mask_fraction = 0)
cat_freq = get_cat_index_freq(d_gen$x)
freq = cat_freq$freq
mu_est = get_cat_mu(cat_freq$freq, n_MC = 5000, old = FALSE)
max(mu_set$est_precis)
3.578706e-08
"


