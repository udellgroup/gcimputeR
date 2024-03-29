% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_core.R
\name{sum_2d_scale}
\alias{sum_2d_scale}
\title{Sum a 2d matrix along one axis}
\usage{
sum_2d_scale(M, v, index = NULL)
}
\arguments{
\item{M}{A 2d matrix with dimension n*k}

\item{index}{a subset of \{1,2,...,n\}}

\item{c}{A vector of length n}
}
\value{
A k*1 matrix
}
\description{
For a 2d matrix M with dimension n*k, sum over dimension 2 with scale c, at entries selected by index in dimension 1
}
\keyword{internal}
