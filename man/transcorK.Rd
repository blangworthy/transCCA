% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcorK.R
\name{transcorK}
\alias{transcorK}
\title{Create transformed Kendall Scatter Matrix}
\usage{
transcorK(x, y, eigenmin = 0.001)
}
\arguments{
\item{x}{The first data matrix to be included in the calculation of the transformed Kendall scatter matrix}

\item{y}{The second data matrix to be included in the calculation of the transformed Kendall scatter matrix}

\item{eigenmin}{The minimum eigenvalue when transforming matrix to be positive definite}
}
\value{
The estimate of the scatter matrix for data using the transformation of Kendall's tau
}
\description{
Creates transformed Kendall's scatter matrix between
two sets of variables, x and y, to be used in scale invariant CCA
}
