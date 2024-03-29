% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/row_function.R
\name{latent_operation_row}
\alias{latent_operation_row}
\title{Operation to conduct in the latent space for each row}
\usage{
latent_operation_row(
  task,
  z,
  lower,
  upper,
  d_index,
  dcat_index,
  corr,
  cat_input = NULL,
  trunc_method = "Iterative",
  n_sample = 5000,
  n_update = 1,
  n_MI = 1
)
}
\arguments{
\item{task}{Task to perform. One of \code{"em", "fillup", "sample"}.}

\item{z}{A row of \code{Z}}

\item{lower}{A row of \code{Lower}}

\item{upper}{A row of \code{Upper}}

\item{d_index}{Boolean vector with \code{TRUE} at ordinal dimensions}

\item{dcat_index}{Boolean vector with \code{TRUE} at categorical dimensions}

\item{corr}{Current copula correlation estimate}

\item{cat_input}{Input for categorical dimensions}

\item{trunc_method}{Method for evaluating truncated normal moments: \code{'Iterative'} or \code{'Sampling'}.}

\item{n_sample}{Number of MC samples, only used when \code{trunc_method} is \code{'Sampling'}}

\item{n_update}{The number of updates, only used when \code{trunc_method} is \code{'Iterative'}}
}
\value{
A list containing
\describe{
  \item{\code{loglik}}{Available when \code{task = 'em'}. The average log-likelihood.}
  \item{\code{Z}}{Available when \code{task = 'em'}. Incomplete \code{Z} with updated observed ordinal entries}
  \item{\code{Zimp}}{Available when \code{task = 'em'} or \code{task == 'fillup'} . Complete \code{Z} with observed entries the same as \code{Z} and missing entries imputed}
  \item{\code{Zimp_sample}}{Available when \code{task = 'sample'}. Multiple imputation samples.}
  \item{\code{C}}{Available when \code{task = 'em'}. The conditional co-variance due to missingness}
  \item{\code{var_ordinal}}{Available when \code{task = 'em'} or \code{task = 'fillup'}. The conditional variance due to truncation, i.e. Var(z|a < z < b)}
}
}
\description{
operation for each latent row
}
\keyword{internal}
