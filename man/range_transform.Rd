% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_fun.R
\name{range_transform}
\alias{range_transform}
\title{Transform original observation X to copula observation Z}
\usage{
range_transform(X, type)
}
\arguments{
\item{X}{Original observation matrix}

\item{type}{Data type. Either \code{"continuous"} or \code{"ordinal"}.}
}
\value{
A list containing
\describe{
  \item{\code{Z}}{Available when \code{type=continuous}. Corresponding value}
  \item{\code{Lower}}{Available when \code{type=ordinal}. Lower boundary of corresponding truncated interval}
  \item{\code{Upper}}{Available when \code{type=ordinal}. Upper boundary of corresponding truncated interval}
}
}
\description{
For continuous \code{X}, compute its corresonding latent value \code{Z}; For ordinal \code{X}, compute its corresponding truncated interval on latent value \code{Z}.
}
\keyword{internal}
