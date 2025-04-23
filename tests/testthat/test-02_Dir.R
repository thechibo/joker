test_that("Dir distr works", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Dir")

  # Errors
  expect_error(Dir(c(-1, 2)))
  expect_error(Dir(-1, 2))

})

test_that("Dir dpqr work", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(r(D)))
  expect_true(is.numeric(d(D, x)))

  # Values
  expect_equal(d(D)(rep(0, 4)), 0)
  expect_equal(sum(x <= 1), 4L * n)
  expect_equal(sum(x >= 0), 4L * n)
  expect_equal(ddir(x[1, ], a, log = TRUE),
               log(ddir(x[1, ], a, log = FALSE)),
               tolerance = 1e-8)


  # 2-Way Calls
  expect_equal(d(D)(1:4 / 10), ddir(1:4 / 10, a))
  expect_equal(d(D)(1:4 / 10), d(D, 1:4 / 10))

  # Error
  expect_error(ddir(x, c(1, 2)))
  expect_error(ddir(x, c(1, 2, 3, -4)))
  expect_error(rdir(-5, 1:4))
  expect_error(rdir(5, c(1, 2, 3, -4)))

})

test_that("Dir moments work", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)

  # Types
  expect_true(is.list(moments(D)))
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(mode(D)))
  expect_true(is.matrix(var(D)))
  expect_true(is.numeric(entro(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Dir likelihood works", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(lldir(x, a)))

  # 2-Way Calls
  expect_equal(lldir(x, a), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

  # ll and lloptim convergence to a0 comparison
  method <- "L-BFGS-B"
  lower <- 1e-5
  upper <- Inf

  par1 <- optim(par = sum(same(D, x)$alpha),
                fn = lloptim,
                gr = dlloptim,
                tx = colMeans(log(x)),
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  par2 <- optim(par = same(D, x)$alpha,
                fn = function(par, x) { ll(Dir(par), x) },
                x = x,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(par1, sum(par2), tolerance = 0.01)

})

test_that("Dir estim works", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(edir(x, type = "mle")))
  expect_true(is.list(edir(x, type = "me")))
  expect_true(is.list(edir(x, type = "same")))

  # 2-Way Calls
  expect_equal(edir(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(edir(x, type = "me"), e(D, x, type = "me"))
  expect_equal(edir(x, type = "same"), e(D, x, type = "same"))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)
  d <- test_consistency("same", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.02)

  # Errors
  expect_error(e(D, x, type = "xxx"))
  expect_error(e(D, x, type = "mle", par0 = "xxx"))

})

test_that("Dir avar works", {

  # Preliminaries
  a <- 1:2
  D1 <- Dir(a)
  D2 <- Beta(1, 2)

  # Types
  expect_true(is.numeric(vdir(a, type = "mle")))
  expect_true(is.numeric(vdir(a, type = "me")))
  expect_true(is.numeric(vdir(a, type = "same")))

  # 2-Way Calls
  expect_equal(vdir(a, type = "mle"), v(D1, type = "mle"))
  expect_equal(vdir(a, type = "me"), v(D1, type = "me"))
  expect_equal(vdir(a, type = "same"), v(D1, type = "same"))
  expect_equal(vdir(a, type = "mle"), avar_mle(D1))
  expect_equal(vdir(a, type = "me"), avar_me(D1))
  expect_equal(vdir(a, type = "same"), avar_same(D1))

  # Dirichlet - Beta comparison
  expect_equal(unname(avar_mle(D1)), unname(avar_mle(D2)), tolerance = 1e-4)
  expect_equal(unname(avar_me(D1)), unname(avar_me(D2)), tolerance = 1e-4)
  expect_equal(unname(avar_same(D1)), unname(avar_same(D2)), tolerance = 1e-4)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Dir small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  a <- 1:4
  D <- Dir(a)
  set.seed(1)

  prm <- list(name = "alpha",
              pos = 1,
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle", "me", "same"),
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

test_that("Dir large metrics work", {

  # Preliminaries
  a <- 1:4
  D <- Dir(a)

  prm <- list(name = "alpha",
              pos = 1,
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- large_metrics(D, prm,
                       est = c("mle", "me", "same"))
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "LargeMetrics")

})
