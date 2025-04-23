test_that("Nbinom distr works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Nbinom")

  # Errors
  expect_error(Nbinom(-10, 0.5))
  expect_error(Nbinom(0, 0.5))
  expect_error(Nbinom(10, 5))
  expect_error(Nbinom(3:4, 0.5))
  expect_error(Nbinom(10, c(0.5, 0.6)))

})

test_that("Nbinom dpqr work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(0), p ^ k, tolerance = 0.01)
  expect_equal(d(D)(1), k * (1 - p) * p ^ k, tolerance = 0.01)
  expect_equal(p(D)(-1), 0)
  expect_equal(qn(D)(1), Inf)
  expect_equal(qn(D)(0), 0)
  expect_equal(sum(x >= 0), n)

  # 2-Way Calls
  expect_equal(d(D)(1), dnbinom(1, k, p))
  expect_equal(p(D)(1), pnbinom(1, k, p))
  expect_equal(qn(D)(1), qnbinom(1, k, p))
  expect_equal(qn(D)(0), qnbinom(0, k, p))
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(1), qn(D, 1))
  expect_equal(qn(D)(0), qn(D, 0))

  # Special Case: Geom
  D1 <- Nbinom(1, 0.7)
  D2 <- Geom(0.7)
  expect_equal(d(D1)(3), d(D2)(3), tolerance = 0.01)
  expect_equal(p(D1)(3), p(D2)(3), tolerance = 0.01)
  expect_equal(qn(D1)(0.7), qn(D2)(0.7), tolerance = 0.01)

})

test_that("Nbinom moments work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)

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

  # Values
  expect_equal(mean(D), k * (1 / p - 1), tolerance = 0.01)

})

test_that("Nbinom likelihood works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llnbinom(x, size = k, prob = p)))

  # 2-Way Calls
  expect_equal(llnbinom(x, k, p), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

})

test_that("Nbinom estim works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(enbinom(x, k, type = "mle")))
  expect_true(is.list(enbinom(x, k, type = "me")))

  # 2-Way Calls
  expect_equal(enbinom(x, k, type = "mle"), e(D, x, type = "mle"))
  expect_equal(enbinom(x, k, type = "me"), e(D, x, type = "me"))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

  # Errors
  expect_error(e(D, x, type = "xxx"))

})

test_that("Nbinom avar works", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vnbinom(k, p, type = "mle")))
  expect_true(is.numeric(vnbinom(k, p, type = "me")))

  # 2-Way Calls
  expect_equal(vnbinom(k, p, type = "mle"), v(D, type = "mle"))
  expect_equal(vnbinom(k, p, type = "me"), v(D, type = "me"))
  expect_equal(vnbinom(k, p, type = "mle"), avar_mle(D))
  expect_equal(vnbinom(k, p, type = "me"), avar_me(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.1)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Nbinom small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

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

test_that("Nbinom large metrics work", {

  # Preliminaries
  k <- 3
  p <- 0.7
  D <- Nbinom(k, p)
  set.seed(1)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

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
