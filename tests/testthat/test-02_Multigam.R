test_that("Multigam distr works", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Multigam")

  # Errors
  expect_error(Multigam(-1, 2))
  expect_error(Multigam(1, -2))

})

test_that("Multigam dpqr work", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(r(D)))
  expect_true(is.numeric(d(D, x)))
  expect_true(is.numeric(dmultigam(x, a, b, log = TRUE)))

  # Values
  expect_equal(d(D)(rep(0, length(a))), 0)
  expect_equal(sum(x < 0), 0L)

  # 2-Way Calls
  expect_equal(d(D)(x[1, ]), dmultigam(x[1, ], shape = a, scale = b))
  expect_equal(d(D)(x[1, ]), d(D, x[1, ]))

  # Errors
  expect_error(dmultigam(x, 1:3, -3))
  expect_error(dmultigam(x, c(1, 2, -3), 3))
  expect_error(dmultigam(x, 1:5, 3))

})

test_that("Multigam moments work", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)

  # Types
  expect_true(is.list(moments(D)))
  expect_true(is.numeric(mean(D)))
  expect_true(is.numeric(var(D)))
  expect_true(is.numeric(finf(D)))

})

test_that("Multigam likelihood works", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llmultigam(x, shape = a, scale = b)))

  # 2-Way Calls
  expect_equal(llmultigam(x, shape = a, scale = b), ll(D, x))
  expect_equal(ll(D)(x), ll(D, x))

  # ll and lloptim convergence to a0 comparison
  method <- "L-BFGS-B"
  lower <- 1e-5
  upper <- Inf

  k <- ncol(x)
  logz <- colMeans(log(fd(x)))
  xk <- mean(x[, k])
  tx <- c(logz, xk)

  par1 <- optim(par = sum(same(D, x)$shape),
               fn = lloptim,
               gr = dlloptim,
               tx = tx,
               distr = D,
               method = method,
               lower = lower,
               upper = upper,
               control = list(fnscale = -1))$par

  b <- xk / par1
  a <- idigamma(logz - log(b))

  par1 <- c(a, b)

  par2 <- optim(par = unlist(same(D, x)),
                fn = function(par, x, distr) {
                  ll(Multigam(par[seq_along(a)], par[length(a) + 1]), x)
                },
                x = x,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  expect_equal(par1, unname(par2), tolerance = 0.01)

})

test_that("Multigam estim works", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(emultigam(x, type = "mle")))
  expect_true(is.list(emultigam(x, type = "me")))
  expect_true(is.list(emultigam(x, type = "same")))

  # 2-Way Calls
  expect_equal(emultigam(x, type = "mle"), e(D, x, type = "mle"))
  expect_equal(emultigam(x, type = "me"), e(D, x, type = "me"))
  expect_equal(emultigam(x, type = "same"), e(D, x, type = "same"))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)
  d <- test_consistency("same", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)

  # Errors
  expect_error(e(D, x, type = "xxx"))
  expect_error(e(D, x, type = "mle", par0 = "xxx"))

})

test_that("Multigam avar works", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)

  # Types
  expect_true(is.numeric(vmultigam(a, b, type = "mle")))
  expect_true(is.numeric(vmultigam(a, b, type = "me")))
  expect_true(is.numeric(vmultigam(a, b, type = "same")))

  # 2-Way Calls
  expect_equal(vmultigam(a, b, type = "mle"), v(D, type = "mle"))
  expect_equal(vmultigam(a, b, type = "me"), v(D, type = "me"))
  expect_equal(vmultigam(a, b, type = "same"), v(D, type = "same"))
  expect_equal(vmultigam(a, b, type = "mle"), avar_mle(D))
  expect_equal(vmultigam(a, b, type = "me"), avar_me(D))
  expect_equal(vmultigam(a, b, type = "same"), avar_same(D))

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Simulations
  d <- test_avar("mle", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.07)
  d <- test_avar("me", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.05)
  d <- test_avar("same", D)
  expect_equal(d$avar_true, d$avar_est, tolerance = 0.07)

  # Errors
  expect_error(v(D, type = "xxx"))

})

test_that("Multigam small metrics work", {

  skip_if(Sys.getenv("JOKER_EXTENDED_TESTS") != "true",
          "Skipping extended test unless JOKER_EXTENDED_TESTS='true'")

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)
  set.seed(1)

  prm <- list(name = "shape",
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

test_that("Multigam large metrics work", {

  # Preliminaries
  a <- 1:3
  b <- 3
  D <- Multigam(a, b)
  set.seed(1)

  prm <- list(name = "shape",
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

  prm <- list(name = "scale",
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
