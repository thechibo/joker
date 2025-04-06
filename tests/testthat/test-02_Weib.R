test_that("Weib distr works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)

  # Types
  expect_s4_class(D, "Distribution")
  expect_s4_class(D, "Weib")

  # Errors
  expect_error(Weib(c(0.1, 0.2, 0.3)))
  expect_error(Weib(-1, 2))
  expect_error(Weib(1, -2))

})

test_that("Weib dpqr work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.function(d(D)))
  expect_true(is.function(p(D)))
  expect_true(is.function(qn(D)))
  expect_true(is.function(r(D)))

  # Values
  expect_identical(d(D)(0), 0)
  expect_identical(p(D)(Inf), 1)
  expect_identical(p(D)(0), 0)
  expect_identical(qn(D)(1), Inf)
  expect_identical(qn(D)(0), 0)
  expect_identical(sum(x >= 0), n)

  # 2-Way Calls
  expect_identical(d(D)(0.4), dweibull(0.4, shape = a, scale = b))
  expect_identical(p(D)(0.4), pweibull(0.4, shape = a, scale = b))
  expect_equal(qn(D)(0.4), qweibull(0.4, shape = a, scale = b),
               tolerance = 1e-8)
  expect_identical(d(D)(0.4), d(D, 0.4))
  expect_identical(p(D)(0.4), p(D, 0.4))
  expect_identical(qn(D)(0.4), qn(D, 0.4))

})

test_that("Weib moments work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)

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

  # Values
  expect_identical(mode(Weib(0.6, 5)), 0)

})

test_that("Weib likelihood works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)
  set.seed(1)
  n <- 1000L
  x <- r(D)(n)

  # Types
  expect_true(is.numeric(llweibull(x, shape = a, scale = b)))

  # 2-Way Calls
  expect_identical(llweibull(x, shape = a, scale = b), ll(D, x))
  expect_identical(ll(D)(x), ll(D, x))

  # ll and lloptim convergence to a0 comparison
  method <- "L-BFGS-B"
  lower <- 1e-5
  upper <- Inf

  par1 <- optim(par = eweibull(x, type = "lme")$shape,
                fn = lloptim,
                gr = dlloptim,
                tx = x,
                distr = D,
                method = method,
                lower = lower,
                upper = upper,
                control = list(fnscale = -1))$par

  par1 <- c(shape = par1, scale = mean(x ^ par1) ^ (1 / par1))

  par2 <- optim(par = unlist(eweibull(x, type = "lme")),
                fn = function(par, x) { ll(Weib(par[1], par[2]), x) },
                x = x,
                method = method,
                lower = c(1e-5, 1e-5),
                upper = c(Inf, Inf),
                control = list(fnscale = -1))$par

  expect_equal(par1, par2, tolerance = 0.01)

})

test_that("Weib estim works", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)
  set.seed(1)
  n <- 100L
  x <- r(D)(n)

  # Types
  expect_true(is.list(eweibull(x, type = "mle")))
  expect_true(is.list(eweibull(x, type = "me")))

  # 2-Way Calls
  expect_identical(eweibull(x, type = "mle"), e(D, x, type = "mle"))
  expect_identical(eweibull(x, type = "me"), e(D, x, type = "me"))

  # Simulations
  d <- test_consistency("me", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)
  d <- test_consistency("mle", D)
  expect_equal(d$prm_true, d$prm_est, tolerance = 0.05)

  # Errors
  expect_error(e(D, x, type = "xxx"))
  expect_error(e(D, x, type = "mle", par0 = "xxx"))
  expect_error(e(D, x, type = "me", par0 = "xxx"))
  expect_warning(e(D, x, type = "me", lower = 0.01))

})

test_that("Weib small metrics work", {

  # Preliminaries
  a <- 2
  b <- 3
  D <- Weib(a, b)
  set.seed(1)

  prm <- list(name = "shape",
              val = seq(0.5, 5, by = 0.5))

  expect_no_error(
    x <- small_metrics(D, prm,
                       est = c("mle"),
                       obs = c(20, 50),
                       sam = 1e2,
                       seed = 1)
  )

  expect_no_error(
    plot(x, save = TRUE, path = tempdir())
  )

  # Types
  expect_s4_class(x, "SmallMetrics")

})
