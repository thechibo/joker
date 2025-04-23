test_that("check_data works", {

  expect_true(is_numatvec(check_data(r(Beta(2, 3), 20))))
  expect_true(is_nummat(check_data(r(Dir(2:4), 20))))
  expect_true(length(check_data(c(1, 2, NA), na.rm = TRUE)) == 2)

  expect_error(check_data(1:10, na.rm = c(TRUE, FALSE)))
  expect_error(check_data(list(1:10)))
  expect_error(check_data(c(1, 2, Inf)))
  expect_error(check_data(c(1, 2, NaN)))
  expect_warning(check_data(c(1, 2, NA), na.rm = FALSE))
  expect_error(check_data(matrix(c(1, 2, NA, 4), 2, 2), na.rm = TRUE))

})

test_that("check_optim works", {

  expect_true(is.character(check_optim("me",
                                       "L-BFGS-B",
                                       1e-5,
                                       Inf,
                                       choices = c("me", "same"),
                                       len = 1)))
  expect_error(check_optim("me",
                           "L-BFGS-B",
                           c(1e-5, 1e-5),
                           c(Inf, Inf),
                           choices = c("me", "same"),
                           len = 1))

  expect_error(check_optim(list("me", "same"),
                           "L-BFGS-B",
                           1e-5,
                           Inf,
                           choices = c("me", "same"),
                           len = 1))

  expect_error(check_optim(10,
                           "L-BFGS-B",
                           "0",
                           Inf,
                           choices = c("me", "same"),
                           len = 1))

  expect_error(check_optim(10,
                           "L-BFGS-B",
                           1e-5,
                           5,
                           choices = c("me", "same"),
                           len = 1))

  expect_error(check_optim(c(1, 2),
                           "L-BFGS-B",
                           1e-5,
                           Inf,
                           choices = c("me", "same"),
                           len = 1))

})

test_that("get functions work", {

  D <- Beta(2, 3)

  expect_true(is.character(get_class_abbr(Beta())))
  expect_identical(get_class_abbr(Beta()), "beta")
  expect_identical(get_class_abbr(Gam()), "gamma")
  expect_identical(get_class_abbr(Stud()), "t")
  expect_identical(get_class_abbr(Fisher()), "f")
  expect_identical(get_class_abbr(Weib()), "weibull")

  expect_s4_class(get_distr_class("bEtA"), "Beta")
  expect_s4_class(get_distr_class("gamma"), "Gam")

  expect_true(is.list(get_params(D, list = TRUE)))
  expect_true(is.atomic(get_params(D, list = FALSE)))

})

test_that("setxofy works", {
  x <- array(1:24, dim = c(2, 3, 4))
  expect_no_error(set2of3(x, 1))
  expect_error(set2of3(x, 6))
})

test_that("seqcol works", {
  expect_error(seqcol(1:5))
  expect_error(seqrow(1:5))
  expect_no_error(seqcol(diag(4)))
  expect_no_error(seqrow(diag(4)))
})

test_that("stats tools work", {

  x <- matrix(1:10, nrow = 2, ncol = 5)

  expect_error(bvar(diag(3)))
  expect_error(bsd(diag(3)))

  expect_error(colVar(1:5))
  expect_error(rowVar(1:5))
  expect_true(is.numeric(rowVar(x)))
  expect_true(is.atomic(rowVar(x)))

})

test_that("inv2x2 works", {
  expect_error(inv2x2(diag(3)))
  expect_error(inv2x2("hello"))
  expect_error(inv2x2(matrix(0, 2, 2)))
  expect_true(is.matrix(inv2x2(diag(2))))
  expect_type(inv2x2(diag(2)), "double")
})

test_that("isit functions work", {

  expect_true(is_atvec(1:5))
  expect_false(is_atvec(list()))

  expect_error(is_symmetric(1:3))
  expect_error(is_symmetric("hello"))
  expect_false(is_symmetric(matrix(1:4, nrow = 2, ncol = 2)))
  expect_true(is_symmetric(diag(3)))

  expect_error(is_pd(1:3))
  expect_error(is_pd("hello"))
  expect_false(is_pd(matrix(1:4, nrow = 2, ncol = 2)))
  expect_false(is_pd(matrix(c(0, 1, -1, 0), 2, 2)))
  expect_false(is_pd(-diag(2)))
  expect_true(is_pd(diag(3)))

  expect_true(is_pos(1:4))
  expect_false(is_pos(c(-1, 2)))
  expect_false(is_pos(c(NaN, 2)))
  expect_false(is_pos(c(Inf, 2)))
  expect_false(is_pos(c(NA, 2)))
})

test_that("fd works", {
  expect_error(fd(list(1:3)))
})
