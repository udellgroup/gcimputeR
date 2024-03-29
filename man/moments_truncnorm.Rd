% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncated_normal.R
\name{moments_truncnorm}
\alias{moments_truncnorm}
\title{Truncated normal mean and variance}
\usage{
moments_truncnorm(mu, std, a, b, tol = 1e-06, mean_only = FALSE)
}
\arguments{
\item{mu}{Mean of the normal distribution}

\item{std}{Standard deviation of the normal distribution}

\item{a}{Left endpoint of the truncated interval}

\item{b}{Right endpoint of the truncated interval}

\item{tol}{Computation is only conducted if the truncated normal probability is larger than \code{tol}}

\item{mean_only}{Only compute the truncated normal variance when `mean_only=FALSE`}
}
\value{
A list containing
\describe{
  \item{\code{mean}}{The truncated normal mean}
  \item{\code{var}}{The truncated normal varaince or `NULL` }
}
}
\description{
Compute the mean and variance of one dimensional truncated normal random variable
}
\keyword{internal}
