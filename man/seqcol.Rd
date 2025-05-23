% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_tools.R
\name{seqcol}
\alias{seqcol}
\alias{seqrow}
\title{Sequence Functions for Matrices}
\usage{
seqcol(x)

seqrow(x)
}
\arguments{
\item{x}{matrix.}
}
\value{
A sequence of integers from 1 to the number of rows or columns of
the matrix.
}
\description{
This set of functions extent the \code{seq_along()} functions for \code{matrix}
objects.
}
\keyword{internal}
