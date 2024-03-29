% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_fun.R
\name{observed_to_latent}
\alias{observed_to_latent}
\title{From observed X to latent Z}
\usage{
observed_to_latent(X, d_index, method = "univariate_mean")
}
\arguments{
\item{X}{Original observation matrix}

\item{d_index}{Boolean vector with \code{TRUE} at ordinal dimensions}

\item{method}{method for initializing the mean}
}
\value{
A list containing
\describe{
  \item{\code{Z}}{Transformed latent matrix}
  \item{\code{Lower}}{Lower boundary for ordinal dimensions. \code{NA} at missing locations.}
  \item{\code{Upper}}{Upper boundary for ordinal dimensions. \code{NA} at missing locations.}
}
}
\description{
Given data matrix \code{X}, prepare the latent matrix \code{Z} and the value bounds for the ordinal dimensions
}
\keyword{internal}
