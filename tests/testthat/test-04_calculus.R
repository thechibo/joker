test_that("Normal calculus works", {

  X <- Norm(0, 1)
  Y <- Norm(4, 2)

  expect_s4_class(X + Y, "Norm")
  expect_s4_class(3 + Y, "Norm")
  expect_s4_class(Y + 3, "Norm")

  expect_s4_class(X - Y, "Norm")
  expect_s4_class(3 - Y, "Norm")
  expect_s4_class(Y - 3, "Norm")

  expect_s4_class(3 * Y, "Norm")
  expect_s4_class(Y * 3, "Norm")
  expect_s4_class(Y / 3, "Norm")

  expect_s4_class(sum(X, Y), "Norm")
  expect_s4_class(exp(Y), "Lnorm")

})
