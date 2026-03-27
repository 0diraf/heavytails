test_that("rpareto returns a numeric vector of the correct length", {
  x <- rpareto(500, alpha = 2, xm = 1)
  expect_type(x, "double")
  expect_length(x, 500)
  expect_true(all(x >= 1))
})

test_that("rpareto returns empty vector for n = 0", {
  x <- rpareto(0, alpha = 2, xm = 1)
  expect_length(x, 0)
})

test_that("pareto_cdf returns values in [0, 1] and matches known points", {
  # F(xmin) = 0, F(2*xmin) = 1 - (1/2)^alpha
  out <- pareto_cdf(c(0.5, 1, 2, 5), xmin = 1, alpha = 2)
  expect_true(all(out >= 0 & out <= 1))
  expect_equal(out[1], 0)                          # x < xmin -> 0
  expect_equal(out[2], 0)                          # x == xmin -> 0
  expect_equal(out[3], 1 - (2/1)^(-2), tolerance = 1e-10)
})

test_that("dpareto returns non-negative values and is zero below xm", {
  out <- dpareto(c(0.5, 1, 2, 10), alpha = 2, xm = 1)
  expect_true(all(out >= 0))
  expect_equal(out[1], 0)   # x < xm
  # density at xm should equal alpha / xm
  expect_equal(out[2], 2 / 1, tolerance = 1e-10)
})

test_that("rpareto validates inputs", {
  expect_error(rpareto(-1, alpha = 2, xm = 1))
  expect_error(rpareto(10, alpha = -1, xm = 1))
  expect_error(rpareto(10, alpha = 2, xm = -1))
  expect_error(rpareto(10, alpha = 2, xm = 0))
})

test_that("pareto_cdf validates inputs", {
  expect_error(pareto_cdf(c(1, 2), xmin = -1, alpha = 2))
  expect_error(pareto_cdf(c(1, 2), xmin = 1, alpha = -1))
  expect_error(pareto_cdf("a", xmin = 1, alpha = 2))
})

test_that("dpareto validates inputs", {
  expect_error(dpareto(c(1, 2), alpha = -1, xm = 1))
  expect_error(dpareto(c(1, 2), alpha = 2, xm = 0))
  expect_error(dpareto("a", alpha = 2, xm = 1))
})
