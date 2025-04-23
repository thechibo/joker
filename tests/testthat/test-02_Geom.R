test_that("Geom distr works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Geom")

  # Errors
  expect_error(Geom(1:2))
  expect_error(Geom(-1))
  expect_error(Geom(4))

})

test_that("Geom dpqr work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_equal(d(D)(-1), 0)
  expect_equal(d(D)(0), p)
  expect_equal(p(D)(0), p)
  expect_equal(p(D)(Inf), 1)
  expect_equal(qn(D)(1), Inf)
  expect_equal(qn(D)(0), 0)
  expect_equal(sum(r(D)(n) >= 0), n)

  # 2-Way Calls
  expect_equal(d(D)(1), dgeom(1, p))
  expect_equal(p(D)(1), pgeom(1, p))
  expect_equal(qn(D)(0.5), qgeom(0.5, p), tolerance = 0.01)
  expect_equal(d(D)(1), d(D, 1))
  expect_equal(p(D)(1), p(D, 1))
  expect_equal(qn(D)(0.5), qn(D, 0.5), tolerance = 0.01)

})

test_that("Geom moments work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

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

  # Warnings
  expect_warning(median(Geom(1 - sqrt(2) / 2)))

})

test_that("Geom likelihood works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llgeom(x, p)))

  # 2-Way Calls
  expect_equal(llgeom(x, p), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

})

test_that("Geom estim works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(eexp(x, type = "mle")))
  expect_true(is.list(eexp(x, type = "me")))

  # 2-Way Calls
  expect_equal(egeom(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(egeom(x, type = "me"), e(D, x, type = "me"))

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

test_that("Geom avar works", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

  # Types
  expect_true(is.numeric(vgeom(p, type = "mle")))
  expect_true(is.numeric(vgeom(p, type = "me")))

  # 2-Way Calls
  expect_equal(vgeom(p, type = "mle"), v(D, type = "mle"))
  expect_equal(vgeom(p, type = "me"), v(D, type = "me"))
  expect_equal(vgeom(p, type = "mle"), avar_mle(D))
  expect_equal(vgeom(p, type = "me"), avar_me(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Geom small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  p <- 0.4
  D <- Geom(p)
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

test_that("Geom large metrics work", {

  # Preliminaries
  p <- 0.4
  D <- Geom(p)

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
