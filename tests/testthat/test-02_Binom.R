test_that("Binom distr works", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Binom")

  # Errors
  expect_error(Binom(-10, 0.5))
  expect_error(Binom(0, 0.5))
  expect_error(Binom(10, 5))
  expect_error(Binom(3:4, 0.5))
  expect_error(Binom(10, c(0.5, 0.6)))

})

test_that("Binom dpqr work", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(N), p ^ N, tolerance = 0.01)
  expect_equal(d(D)(0), (1 - p) ^ N, tolerance = 0.01)
  expect_equal(p(D)(N), 1)
  expect_equal(qn(D)(1), N)
  expect_equal(qn(D)(0), 0)
  expect_equal(sum(x %in% 0:N), n)

  # 2-Way Calls
  expect_equal(d(D)(1), dbinom(1, N, p))
  expect_equal(p(D)(1), pbinom(1, N, p))
  expect_equal(qn(D)(1), qbinom(1, N, p))
  expect_equal(qn(D)(0), qbinom(0, N, p))
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(1), qn(D, 1))
  expect_equal(qn(D)(0), qn(D, 0))

})

test_that("Binom moments work", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(sd(D)))
  expect_true(is.numeric(skew(D)))
  expect_true(is.numeric(kurt(D)))
  expect_true(is.numeric(finf(D)))

  # Warnings
  expect_warning(moments(D))
  expect_warning(entro(D))

  # Values
  expect_equal(mean(D), N * p)
  expect_equal(var(D), N * p * (1 - p))

})

test_that("Binom likelihood works", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llbinom(x, size = N, prob = p)))

  # 2-Way Calls
  expect_equal(llbinom(x, N, p), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

})

test_that("Binom estim works", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(ebinom(x, size = N, type = "mle")))
  expect_true(is.list(ebinom(x, size = N, type = "me")))

  # 2-Way Calls
  expect_equal(ebinom(x, N, type = "mle"), e(D, x, type = "mle"),
                   tolerance = 1e-16)
  expect_equal(ebinom(x, N, type = "me"), e(D, x, type = "me"),
                   tolerance = 1e-16)

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

  # Errors
  expect_error(mle(Binom(1, 0.5), c(3, 5, 4)))
  expect_error(e(D, x, type = "xxx"))

})

test_that("Binom avar works", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vbinom(N, p, type = "mle")))
  expect_true(is.numeric(vbinom(N, p, type = "me")))

  # 2-Way Calls
  expect_equal(vbinom(N, p, type = "mle"), v(D, type = "mle"),
                   tolerance = 1e-16)
  expect_equal(vbinom(N, p, type = "me"), v(D, type = "me"),
                   tolerance = 1e-16)
  expect_equal(vbinom(N, p, type = "mle"), avar_mle(D))
  expect_equal(vbinom(N, p, type = "me"), avar_me(D))

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

test_that("Binom small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
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

test_that("Binom large metrics work", {

  # Preliminaries
  N <- 10
  p <- 0.7
  D <- Binom(N, p)
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
