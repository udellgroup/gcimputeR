% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_fun.R
\name{get_cat_bounds}
\alias{get_cat_bounds}
\title{Get the lower and upper bounds for truncated normal moments calculation at categorical columns}
\usage{
get_cat_bounds(X_cat, mu, cat_index_list, check = FALSE, old = FALSE)
}
\arguments{
\item{X_cat}{(incomplete) categorical data matrix}

\item{mu}{mean vector for categorical}

\item{cat_index_list}{list indicating the number of levels for each categorical variable}

\item{check}{Sanity check on the output?}

\item{old}{Use previous stringent categorical model?}
}
\value{
Imputed data matrix
}
\description{
Impute orginal incomplete observation using imputed \code{Z} matrix and marginal distribution information drawn from original observation \code{X}
}
\keyword{internal}
