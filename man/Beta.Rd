% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Beta.R
\name{Beta}
\alias{Beta}
\alias{d,Beta,numeric-method}
\alias{p,Beta,numeric-method}
\alias{qn,Beta,numeric-method}
\alias{r,Beta,numeric-method}
\alias{mean,Beta-method}
\alias{median,Beta-method}
\alias{mode,Beta-method}
\alias{var,Beta-method}
\alias{sd,Beta-method}
\alias{skew,Beta-method}
\alias{kurt,Beta-method}
\alias{entro,Beta-method}
\alias{finf,Beta-method}
\alias{llbeta}
\alias{ll,Beta,numeric-method}
\alias{ebeta}
\alias{mle,Beta,numeric-method}
\alias{me,Beta,numeric-method}
\alias{same,Beta,numeric-method}
\alias{vbeta}
\alias{avar_mle,Beta-method}
\alias{avar_me,Beta-method}
\alias{avar_same,Beta-method}
\title{Beta Distribution}
\usage{
Beta(shape1 = 1, shape2 = 1)

\S4method{d}{Beta,numeric}(distr, x, log = FALSE)

\S4method{p}{Beta,numeric}(distr, q, lower.tail = TRUE, log.p = FALSE)

\S4method{qn}{Beta,numeric}(distr, p, lower.tail = TRUE, log.p = FALSE)

\S4method{r}{Beta,numeric}(distr, n)

\S4method{mean}{Beta}(x)

\S4method{median}{Beta}(x)

\S4method{mode}{Beta}(x)

\S4method{var}{Beta}(x)

\S4method{sd}{Beta}(x)

\S4method{skew}{Beta}(x)

\S4method{kurt}{Beta}(x)

\S4method{entro}{Beta}(x)

\S4method{finf}{Beta}(x)

llbeta(x, shape1, shape2)

\S4method{ll}{Beta,numeric}(distr, x)

ebeta(x, type = "mle", ...)

\S4method{mle}{Beta,numeric}(
  distr,
  x,
  par0 = "same",
  method = "L-BFGS-B",
  lower = 1e-05,
  upper = Inf,
  na.rm = FALSE
)

\S4method{me}{Beta,numeric}(distr, x, na.rm = FALSE)

\S4method{same}{Beta,numeric}(distr, x, na.rm = FALSE)

vbeta(shape1, shape2, type = "mle")

\S4method{avar_mle}{Beta}(distr)

\S4method{avar_me}{Beta}(distr)

\S4method{avar_same}{Beta}(distr)
}
\arguments{
\item{shape1, shape2}{numeric. The non-negative distribution parameters.}

\item{distr}{an object of class \code{Beta}.}

\item{x}{For the density function, \code{x} is a numeric vector of quantiles. For
the moments functions, \code{x} is an object of class \code{Beta}. For the
log-likelihood and the estimation functions, \code{x} is the sample of
observations.}

\item{log, log.p}{logical. Should the logarithm of the probability be
returned?}

\item{q}{numeric. Vector of quantiles.}

\item{lower.tail}{logical. If TRUE (default), probabilities are
\eqn{P(X \leq x)}, otherwise \eqn{P(X > x)}.}

\item{p}{numeric. Vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to
be the number required.}

\item{type}{character, case ignored. The estimator type (mle, me, or same).}

\item{...}{extra arguments.}

\item{par0, method, lower, upper}{arguments passed to optim for the mle
optimization. See Details.}

\item{na.rm}{logical. Should the \code{NA} values be removed?}
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
The Beta distribution is an absolute continuous probability distribution with
support \eqn{S = [0,1]}, parameterized by two shape parameters,
\eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\details{
The probability density function (PDF) of the Beta distribution is given by:
\deqn{ f(x; \alpha, \beta) = \frac{x^{\alpha - 1} (1 - x)^{\beta -
1}}{B(\alpha, \beta)},
\quad \alpha\in\mathbb{R}_+, \, \beta\in\mathbb{R}_+,}
for \eqn{x \in S = [0, 1]}, where \eqn{B(\alpha, \beta)} is the Beta
function:
\deqn{ B(\alpha, \beta) = \int_0^1 t^{\alpha - 1} (1 - t)^{\beta - 1} dt.}

The MLE of the beta distribution parameters is not available in closed form
and has to be approximated numerically. This is done with \code{optim()}.
Specifically, instead of solving a bivariate optimization problem w.r.t
\eqn{(\alpha, \beta)}, the optimization can be performed on the parameter
sum \eqn{\alpha_0:=\alpha + \beta \in(0,+\infty)}. The default method used
is the L-BFGS-B method with lower bound \code{1e-5} and upper bound \code{Inf}. The
\code{par0} argument can either be a numeric (satisfying \verb{lower <= par0 <= upper})
or a character specifying the closed-form estimator to be used as
initialization for the algorithm (\code{"me"} or \code{"same"} - the default value).
}
\examples{
# -----------------------------------------------------
# Beta Distribution Example
# -----------------------------------------------------

# Create the distribution
a <- 3
b <- 5
D <- Beta(a, b)

# ------------------
# dpqr Functions
# ------------------

d(D, c(0.3, 0.8, 0.5)) # density function
p(D, c(0.3, 0.8, 0.5)) # distribution function
qn(D, c(0.4, 0.8)) # inverse distribution function
x <- r(D, 100) # random generator function

# alternative way to use the function
df <- d(D) ; df(x) # df is a function itself

# ------------------
# Moments
# ------------------

mean(D) # Expectation
var(D) # Variance
sd(D) # Standard Deviation
skew(D) # Skewness
kurt(D) # Excess Kurtosis
entro(D) # Entropy
finf(D) # Fisher Information Matrix

# List of all available moments
mom <- moments(D)
mom$mean # expectation

# ------------------
# Point Estimation
# ------------------

ll(D, x)
llbeta(x, a, b)

ebeta(x, type = "mle")
ebeta(x, type = "me")
ebeta(x, type = "same")

mle(D, x)
me(D, x)
same(D, x)
e(D, x, type = "mle")

mle("beta", x) # the distr argument can be a character

# ------------------
# Estimator Variance
# ------------------

vbeta(a, b, type = "mle")
vbeta(a, b, type = "me")
vbeta(a, b, type = "same")

avar_mle(D)
avar_me(D)
avar_same(D)

v(D, type = "mle")
}
\references{
\itemize{
\item Tamae, H., Irie, K. & Kubokawa, T. (2020), A score-adjusted approach to
closed-form estimators for the gamma and beta distributions, Japanese Journal
of Statistics and Data Science 3, 543–561.
\item Papadatos, N. (2022), On point estimators for gamma and beta distributions,
arXiv preprint arXiv:2205.10799.
}
}
\seealso{
Functions from the \code{stats} package: \code{\link[=dbeta]{dbeta()}}, \code{\link[=pbeta]{pbeta()}}, \code{\link[=qbeta]{qbeta()}},
\code{\link[=rbeta]{rbeta()}}
}
