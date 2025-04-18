% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Fisher.R
\name{Fisher}
\alias{Fisher}
\alias{d,Fisher,numeric-method}
\alias{p,Fisher,numeric-method}
\alias{qn,Fisher,numeric-method}
\alias{r,Fisher,numeric-method}
\alias{mean,Fisher-method}
\alias{median,Fisher-method}
\alias{mode,Fisher-method}
\alias{var,Fisher-method}
\alias{sd,Fisher-method}
\alias{skew,Fisher-method}
\alias{kurt,Fisher-method}
\alias{entro,Fisher-method}
\alias{llf}
\alias{ll,Fisher,numeric-method}
\title{Fisher Distribution}
\usage{
Fisher(df1 = 1, df2 = 1)

\S4method{d}{Fisher,numeric}(distr, x, log = FALSE)

\S4method{p}{Fisher,numeric}(distr, q, lower.tail = TRUE, log.p = FALSE)

\S4method{qn}{Fisher,numeric}(distr, p, lower.tail = TRUE, log.p = FALSE)

\S4method{r}{Fisher,numeric}(distr, n)

\S4method{mean}{Fisher}(x)

\S4method{median}{Fisher}(x)

\S4method{mode}{Fisher}(x)

\S4method{var}{Fisher}(x)

\S4method{sd}{Fisher}(x)

\S4method{skew}{Fisher}(x)

\S4method{kurt}{Fisher}(x)

\S4method{entro}{Fisher}(x)

llf(x, df1, df2)

\S4method{ll}{Fisher,numeric}(distr, x)
}
\arguments{
\item{df1, df2}{numeric. The distribution degrees of freedom parameters.}

\item{distr}{an object of class \code{Fisher}.}

\item{x}{For the density function, \code{x} is a numeric vector of quantiles. For
the moments functions, \code{x} is an object of class \code{Fisher}. For the
log-likelihood functions, \code{x} is the sample of observations.}

\item{log, log.p}{logical. Should the logarithm of the probability be
returned?}

\item{q}{numeric. Vector of quantiles.}

\item{lower.tail}{logical. If TRUE (default), probabilities are
\eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.}

\item{p}{numeric. Vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to
be the number required.}
}
\value{
Each type of function returns a different type of object:
\itemize{
\item Distribution Functions: When supplied with one argument (\code{distr}), the
\code{d()}, \code{p()}, \code{q()}, \code{r()}, \code{ll()} functions return the density, cumulative
probability, quantile, random sample generator, and log-likelihood functions,
respectively. When supplied with both arguments (\code{distr} and \code{x}), they
evaluate the aforementioned functions directly.
\item Moments: Returns a numeric, either vector or matrix depending on the moment
and the distribution. The \code{moments()} function returns a list with all the
available methods.
\item Estimation: Returns a list, the estimators of the unknown parameters. Note
that in distribution families like the binomial, multinomial, and negative
binomial, the size is not returned, since it is considered known.
\item Variance: Returns a named matrix. The asymptotic covariance matrix of the
estimator.
}
}
\description{
The Fisher (F) distribution is an absolute continuous probability
distribution that arises frequently in the analysis of variance (ANOVA) and
in hypothesis testing. It is defined by two degrees of freedom parameters
\eqn{d_1 > 0} and \eqn{d_2 > 0}.
}
\details{
The probability density function (PDF) of the F-distribution is given by:
\deqn{ f(x; d_1, d_2) = \frac{\sqrt{\left(\frac{d_1 x}{d_1 x +
d_2}\right)^{d_1} \left(\frac{d_2}{d_1 x + d_2}\right)^{d_2}}}{x B(d_1/2,
d_2/2)}, \quad x > 0 .}
}
\examples{
# -----------------------------------------------------
# Fisher Distribution Example
# -----------------------------------------------------

# Create the distribution
df1 <- 14 ; df2 <- 20
D <- Fisher(df1, df2)

# ------------------
# dpqr Functions
# ------------------

d(D, c(0.3, 2, 10)) # density function
p(D, c(0.3, 2, 10)) # distribution function
qn(D, c(0.4, 0.8)) # inverse distribution function
x <- r(D, 100) # random generator function

# alternative way to use the function
df <- d(D) ; df(x) # df is a function itself

# ------------------
# Moments
# ------------------

mean(D) # Expectation
median(D) # Median
mode(D) # Mode
var(D) # Variance
sd(D) # Standard Deviation
skew(D) # Skewness
kurt(D) # Excess Kurtosis
entro(D) # Entropy

# List of all available moments
mom <- moments(D)
mom$mean # expectation

# ------------------
# Point Estimation
# ------------------

ll(D, x)
llf(x, df1, df2)

}
\seealso{
Functions from the \code{stats} package: \code{\link[=df]{df()}}, \code{\link[=pf]{pf()}}, \code{\link[=qf]{qf()}}, \code{\link[=rf]{rf()}}
}
