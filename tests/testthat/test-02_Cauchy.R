test_that("Cauchy distr works", {

  # Preliminaries
  D <- Cauchy(0.7)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Cauchy")

  # Errors
  expect_error(Cauchy(c(0.1, 0.2)))
  expect_error(Cauchy(2, - 1))

})

test_that("Cauchy dpqr work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(m), 1 / (pi * g))
  expect_equal(p(D)(m), 0.5)
  expect_equal(qn(D)(0.5), m)
  expect_equal(qn(D)(0), -Inf)
  expect_equal(qn(D)(1), Inf)

  # 2-Way Calls
  expect_equal(d(D)(1), dcauchy(1, m, g))
  expect_equal(p(D)(1), pcauchy(1, m, g))
  expect_equal(qn(D)(1), qcauchy(1, m, g))
  expect_equal(qn(D)(0), qcauchy(0, m, g))
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(1), qn(D, 1))
  expect_equal(qn(D)(0), qn(D, 0))

})

test_that("Cauchy moments work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale

  # Types
  expect_warning(mean(D))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_warning(var(D))
  expect_warning(sd(D))
  expect_warning(skew(D))
  expect_warning(kurt(D))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Cauchy likelihood works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llcauchy(x, m, g)))

  # 2-Way Calls
  expect_equal(llcauchy(x, m, g), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

  # ll and lloptim convergence comparison
  method <- "L-BFGS-B"
  lower <- c(-Inf, 1e-5)
  upper <- c(Inf, Inf)

  par1 <- optim(par = me(D, x),
                fn = lloptim,
                gr = dlloptim,
                tx = x,
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  par2 <- optim(par = me(D, x),
                fn = function(par, x) { ll(Cauchy(par[1], par[2]), x) },
                x = x,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(par1, par2, tolerance = 0.01)

})

test_that("Cauchy estim works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(ecauchy(x, type = "mle")))
  expect_true(is.list(ecauchy(x, type = "me")))

  # 2-Way Calls
  expect_equal(ecauchy(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(ecauchy(x, type = "me"), e(D, x, type = "me"))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.03)

  # Errors
  expect_error(e(D, x, type = "xxx"))
  expect_error(e(D, x, type = "mle", par0 = "xxx"))

})

test_that("Cauchy avar works", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  expect_true(is.matrix(vcauchy(m, g)))

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Cauchy small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "location",
              val = seq(-2, 2, by = 1))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1,
                       bar = FALSE)
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "SmallMetrics")

})

## NEEDS CHECKING
test_that("Cauchy large metrics work", {

  # Preliminaries
  D <- Cauchy(2, 1)
  m <- D@location
  g <- D@scale
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "location",
              val = seq(-2, 2, by = 1))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

  prm <- list(name = "scale",
              val = seq(0.5, 2, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

})
