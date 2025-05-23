test_that("Bern distr works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Bern")

  # Errors
  expect_error(Bern(c(0.1, 0.2)))
  expect_error(Bern(2))
  expect_error(Bern(-1))

})

test_that("Bern dpqr work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(1), p)
  expect_equal(p(D)(1), 1)
  expect_equal(qn(D)(1), 1)
  expect_equal(qn(D)(0), 0)
  expect_equal(sum(x %in% c(0, 1)), n)

  # 2-Way Calls
  expect_equal(d(D)(1), dbern(1, p))
  expect_equal(p(D)(1), pbern(1, p))
  expect_equal(qn(D)(1), qbern(1, p))
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(1), qn(D, 1))

})

test_that("Bern moments work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)

  # Types
  expect_true(is.list(moments(D)))
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
  expect_equal(mean(D), p)
  expect_equal(var(D), p * (1 - p))
  expect_equal(median(Bern(0.3)), 0)
  expect_equal(mode(Bern(0.3)), 0)

  # Warnings
  expect_warning(median(Bern(0.5)))
  expect_warning(mode(Bern(0.5)))

})

test_that("Bern likelihood works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llbern(x, p)))

  # 2-Way Calls
  expect_equal(llbern(x, p), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

})

test_that("Bern estim works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(ebern(x, type = "mle")))
  expect_true(is.list(ebern(x, type = "me")))

  # 2-Way Calls
  expect_equal(ebern(x, type = "mle"), e(D, x, type = "mle"),
                   tolerance = 1e-16)
  expect_equal(ebern(x, type = "me"), e(D, x, type = "me"),
                   tolerance = 1e-16)

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.01)

  # Errors
  expect_error(e(D, type = "xxx"))

})

test_that("Bern avar works", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(vbern(p, type = "mle")))
  expect_true(is.numeric(vbern(p, type = "me")))

  # 2-Way Calls
  expect_equal(vbern(p, type = "mle"), v(D, type = "mle"),
                   tolerance = 1e-16)
  expect_equal(vbern(p, type = "me"), v(D, type = "me"))
  expect_equal(vbern(p, type = "mle"), avar_mle(D))
  expect_equal(vbern(p, type = "me"), avar_me(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est["prob"], tolerance = 0.01)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est["prob"], tolerance = 0.01)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Bern small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
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

test_that("Bern large metrics work", {

  # Preliminaries
  p <- 0.7
  D <- Bern(p)
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
