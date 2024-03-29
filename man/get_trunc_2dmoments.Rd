% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncated_normal.R
\name{get_trunc_2dmoments}
\alias{get_trunc_2dmoments}
\title{Compute multivariate truncated normal mean and cov}
\usage{
get_trunc_2dmoments(mean, cov, lower, upper, n_sample = 5000)
}
\arguments{
\item{mean}{Normal mean vector}

\item{cov}{Normal covariance matrix}

\item{lower}{Lower boundary of truncated intervals}

\item{upper}{Upper boundary of truncated intervals}

\item{n_sample}{Number of samples to use for sampling methods}
}
\value{
A list containing
\describe{
  \item{\code{mean}}{Mean vector}
  \item{\code{cov}}{Covariance matrix}
}
}
\description{
compute multivariate truncated normal mean and cov by sampling
}
\keyword{internal}
