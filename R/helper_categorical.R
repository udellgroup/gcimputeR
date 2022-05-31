relabel <- function(x, label){
  stopifnot(check_cat_label(x))
  x = as.factor(x)
  levels(x) = label
  as.numeric(levels(x))[x]
}

cat_index_from_list <- function(l) as.integer(names(l))

nominal_z_to_x_col <- function(z, old = FALSE){
  argmax = which.max(z)
  if (old){
    if (z[argmax]<0) argmax = 1 else argmax = argmax + 1
  }
  argmax
}

#' Encode category names to integers
#'
#' @description For a category variable with k categories, rename the category names to 1,2,...,k in the order of appearance.
#' @param x samples of  a categorical variable
#' @return A list containing
#' \describe{
#'   \item{\code{x}}{Integer encoded \code{x}}
#'   \item{\code{xlevels}}{Original category names corresponding to 1,2,...,k.}
#' }
#' @export
cat_to_integers <- function(x){
  x = droplevels(as.factor(x))
  xlevels = levels(x)
  nlevel = nlevels(x)
  levels(x) = 1:nlevel
  list(x=as.numeric(x), xlevels = xlevels)
}

check_cat_label <- function(x){
  x = x[!is.na(x)]
  xmax = max(x)
  xmin = min(x)
  nlevel = nlevels(as.factor(x))
  xmin == 1 & xmax == nlevel
}

adjust_index_list <- function(index_list){
  start = 1
  for (i in seq_along(index_list)){
    vals = index_list[[i]]
    index_list[[i]] = vals - vals[1] + start
    start = start + length(vals)
  }
  index_list
}

get_cat_slicing_index <- function(x_cat, cat_index_list, keep = 'observed', d_cat=NULL){
  if (is.character(keep)){
    if (keep == 'observed') index_incat = !is.na(x_cat)
    else if (keep == 'missing') index_incat = is.na(x_cat)
    else stop('invalid char keep')
  }else if (is.logical(keep) & length(keep) == length(x_cat)){
    index_incat = keep
  } else stop('invalid keep')

  if (is.null(d_cat)) d_cat= sum(purrr::map_int(cat_index_list, length))
  index_cat = logical(d_cat)
  if (any(index_incat)){
    intindex_cat = purrr::reduce(cat_index_list[index_incat], c)
    index_cat[intindex_cat] = TRUE
  }

  list('incat'=index_incat, 'cat'=index_cat)
}



create_cat_index_list <- function(cat_index_level){
  cat_index_list = list()
  start = 1
  for (i in seq_along(cat_index_level)){
    l = cat_index_level[i]
    if (!is.na(l)){
      cat_index_list[[as.character(i)]] = start:(start+l-1)
      start = start+l
    }else start = start+1
  }
  cat_index_list
}


#' Copula correlation projection
#'
#' @description Project a copula correlation matrix to be the identity at each categorical block
#' @param sigma copula  correlation
#' @param eps minimal allowed eigenvalue
#' @inheritParams  get_cat_bounds
#' @export
project_to_nominal_corr <- function(sigma, cat_index_list, eps=1e-5){
  p = ncol(sigma)
  A = diag(nrow = p, ncol = p)
  for (i in seq_along(cat_index_list)){
    index = cat_index_list[[i]]
    o_eigen = eigen(sigma[index,index])
    #m = diag(1/sqrt(o_eigen$values)) %*% t(o_eigen$vectors)
    eigen_values = o_eigen$values
    if (min(eigen_values) < 1e-5){
      warning('Projection skipped: small eigenvalue in a categorical block')
      A[index,index] = diag(length(index))
    }else{
      m = o_eigen$vectors %*% sqrt(diag(1/eigen_values)) %*% t(o_eigen$vectors)
      A[index,index] = m
    }
  }

  sigma = A %*% sigma %*% t(A)

  for (i in seq_along(cat_index_list)){
    index = cat_index_list[[i]]
    l = length(index)
    sigma[index,index] = diag(nrow = l, ncol = l)
  }
  sigma
}

"
project_to_nominal_corr_simple <- function(sigma, cat_index_list, eps = 1e-3){
  p = ncol(sigma)
  for (i in seq_along(cat_index_list)){
    index = cat_index_list[[i]]
    sigma[index,index] = diag(length(index))
  }
  o_eigen = eigen(sigma)$values
  if (min(o_eigen)<0){
    sigma = sigma - diag(p) * (min(o_eigen)-eps)
    sigma = cov2cor(sigma)
  }
  sigma
}
"



#' Transform latent for cat
#'
#' @description  Transform latent matrix corresponding to latent categorical variables so that after transformation, interval truncated gaussian is achieved
#' @param Z Matrix corresponding to ONLY latent categorical variables
#' @inheritParams get_cat_bounds
#' @keywords internal
Z_to_original_trunc <- function(Z, X_cat, cat_index_list, old = FALSE){
  n = nrow(Z)
  for (i in 1:n){
    z = Z[i,]
    x_cat = X_cat[i,]
    obs_indices = !is.na(z)
    cat_obs = !is.na(x_cat)
    if (any(cat_obs)){
      A = x_to_A(x = x_cat[cat_obs], cat_index_list = cat_index_list[cat_obs], old = old)
      if (!is.null(A)) z[obs_indices] = A %*% z[obs_indices]
      Z[i,] = z
    }
  }
  Z
}


# Find matrix A such that A^2=I and Ax is interval truncated. The algorithm forms A = diag(A1,...,Ap)
# for x=(x1,...,xp).
# x must not have missing entries
x_to_A <- function(x, cat_index_list, d_cat=NULL, adjust=TRUE, test=TRUE, old = FALSE){
  if (adjust) cat_index_list = adjust_index_list(cat_index_list)
  if (test){
    if (!is.null(d_cat)){
      if (d_cat != sum(purrr::map_int(cat_index_list, length))) stop('invalid input')
    }
  }
  if (is.null(d_cat)) d_cat = sum(purrr::map_int(cat_index_list, length))
  if (any(is.na(x))) stop('invalid x')

  if (old){
    index_notbase = get_cat_slicing_index(x, cat_index_list, keep = x!=1, d_cat = d_cat)
    if (any(index_notbase$incat)){
      A = diag(nrow = d_cat, ncol = d_cat)
      # for each xi != 1
      for (i in which(index_notbase$incat)){
        index = cat_index_list[[i]]
        x_index = x[i]-1
        Ai = -diag(length(index))
        Ai[,x_index] = 1
        A[index,index] = Ai
      }
    }else{
      A = NULL
    }
  }else{
    A = diag(nrow = d_cat, ncol = d_cat)
    p = length(x)
    for (i in 1:p){
      index = cat_index_list[[i]]
      Ai = -diag(length(index))
      Ai[,x[i]] = 1
      A[index,index] = Ai
    }
  }
  A
}

# For cov(x)=sigma and x contains xcat as a subvector, compute cov(xnew) with xcat replaced by A*xcat
A_sigma_tA_at_cat <- function(sigma, A, cat_index = NULL, A_at_left = TRUE){
  if (is.null(A)) return(sigma)
  p = ncol(sigma)
  if (is.null(cat_index)){
    A_all = A
  }else{
    if (length(cat_index) != p) stop('invalid cat_index and sigma')
    A_all = diag(p)
    A_all[cat_index,cat_index] = A
  }
  if (A_at_left) A_all %*% sigma %*% t(A_all)
  else  t(A_all) %*% sigma %*% A_all
}
