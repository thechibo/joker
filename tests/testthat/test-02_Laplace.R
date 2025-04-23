test_that("Laplace distr works", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Laplace")

  # Errors
  expect_error(Laplace(c(0, 1), 2))
  expect_error(Laplace(0, -1))

})

test_that("Laplace dpqr work", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))
  expect_true(is.numeric(d(D, x, log = TRUE)))
  expect_true(is.numeric(p(D, x, log.p = TRUE)))
  expect_true(is.numeric(p(D, x, lower.tail = FALSE)))
  expect_true(is.numeric(qn(D, log(0.8), log.p = TRUE)))
  expect_true(is.numeric(qn(D, 0.8, lower.tail = FALSE)))

  # Values
  expect_equal(p(D)(mu), 0.5)
  expect_equal(p(D)(Inf), 1)
  expect_equal(qn(D)(1), Inf)
  expect_equal(qn(D)(0.5), mu)
  expect_equal(qn(D)(0), -Inf)

  # 2-Way Calls
  expect_equal(d(D)(1), dlaplace(1, mu, sigma))
  expect_equal(p(D)(1), plaplace(1, mu, sigma))
  expect_equal(qn(D)(0.5), qlaplace(0.5, mu, sigma), tolerance = 0.01)
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(0.5), qn(D)(0.5), tolerance = 0.01)

})

test_that("Laplace moments work", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(median(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Laplace likelihood works", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(lllaplace(x, mu, sigma)))

  # 2-Way Calls
  expect_equal(lllaplace(x, mu, sigma), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

})

test_that("Laplace estim works", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(elaplace(x, type = "mle")))
  expect_true(is.list(elaplace(x, type = "me")))

  # 2-Way Calls
  expect_equal(elaplace(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(elaplace(x, type = "me"), e(D, x, type = "me"))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

  # Errors
  expect_error(e(D, x, type = "xxx"))

})

test_that("Laplace avar works", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vlaplace(mu, sigma, type = "mle")))
  expect_true(is.numeric(vlaplace(mu, sigma, type = "me")))

  # 2-Way Calls
  expect_equal(vlaplace(mu, sigma, type = "mle"), v(D, type = "mle"))
  expect_equal(vlaplace(mu, sigma, type = "me"), v(D, type = "me"))
  expect_equal(vlaplace(mu, sigma, type = "mle"), avar_mle(D))
  expect_equal(vlaplace(mu, sigma, type = "me"), avar_me(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D, n = 2e4, m = 2e3)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D, n = 2e4, m = 2e3)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Laplace small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "mu",
              val = seq(0.5, 5, by = 0.5))

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

  prm <- list(name = "sigma",
              val = seq(0.5, 5, by = 0.5))

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

test_that("Laplace large metrics work", {

  # Preliminaries
  mu <- 3
  sigma <- 1
  D <- Laplace(mu, sigma)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  prm <- list(name = "mu",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

  prm <- list(name = "sigma",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

})
