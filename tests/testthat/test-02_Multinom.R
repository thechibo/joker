test_that("Multinom distr works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Multinom")

  # Errors
  expect_error(Multinom(-10, 0.5))
  expect_error(Multinom(0, 0.5))
  expect_error(Multinom(10, 5))
  expect_error(Multinom(3:4, 0.5))

})

test_that("Multinom dpqr work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(r(D)))
  expect_true(is.numeric(d(D, x[, 1])))

  # Values
  expect_equal(d(D)(c(N, 0, 0)), p[1] ^ N, tolerance = 0.01)
  expect_equal(sum(x %in% 0:N), length(p)*n)
  expect_equal(sum(colSums(x) == N), n)

  # 2-Way Calls
  expect_equal(d(Multinom(N, p))(c(4, 3, 3)),
                   dmultinom(c(4, 3, 3), N, p))
  expect_equal(d(Multinom(N, p))(c(4, 3, 3)),
                   d(Multinom(N, p), c(4, 3, 3)))

  # Special Case: Binomial
  expect_equal(d(Multinom(N, c(0.3, 0.7)))(c(2, N-2)),
                   dbinom(2, N, 0.3), tolerance = 0.01)

})

test_that("Multinom moments work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  # Types
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))
  expect_true(is.numeric(finf(Multinom(N, c(0.6, 0.4)))))

  # Values
  expect_equal(mean(D), N * p)
  expect_equal(var(D)[1, 1], N * p[1] * (1 - p[1]), tolerance = 0.01)

})

test_that("Multinom likelihood works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llmultinom(x, size = N, prob = p)))

  # 2-Way Calls
  expect_equal(llmultinom(x, N, p), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

  # Error
  x[1, 1] <- x[1, 1] - 1
  expect_error(ll(D, x))

})

test_that("Multinom estim works", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(emultinom(x, type = "mle")))
  expect_true(is.list(emultinom(x, type = "me")))

  # 2-Way Calls
  expect_equal(emultinom(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(emultinom(x, type = "me"), e(D, x, type = "me"))

  # Error
  x[1, 1] <- x[1, 1] - 1
  expect_error(mle(D, x))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

  # Errors
  expect_error(e(D, x, type = "xxx"))

})

test_that("Multinom avar works", {

  # Preliminaries
  N <- 1e3
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)
  k <- length(p)

  # Types
  expect_true(is.numeric(vmultinom(N, p, type = "mle")))
  expect_true(is.numeric(vmultinom(N, p, type = "me")))

  # 2-Way Calls
  expect_equal(vmultinom(N, p, type = "mle"), v(D, type = "mle"))
  expect_equal(vmultinom(N, p, type = "me"), v(D, type = "me"))
  expect_equal(vmultinom(N, p, type = "mle"), avar_mle(D))
  expect_equal(vmultinom(N, p, type = "me"), avar_me(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.1)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est[1:(k-1), 1:(k-1)], tolerance = 0.1)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Multinom small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

  expect_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1,
                       bar = FALSE)
  )

})

test_that("Multinom large metrics work", {

  # Preliminaries
  N <- 10
  p <- c(0.7, 0.2, 0.1)
  D <- Multinom(N, p)

  prm <- list(name = "prob",
              val = seq(0.5, 0.8, by = 0.1))

  expect_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me"))
  )

})
