test_that("lr_test_pareto returns a data.frame with the correct columns", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  out <- lr_test_pareto(x, xm = 1)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("alternative", "ll_pareto", "ll_alternative",
                    "lr_statistic", "p_value", "preferred") %in% names(out)))
  expect_equal(nrow(out), 3)   # exponential, lognormal, weibull
})

test_that("lr_test_pareto prefers pareto on genuine Pareto data", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  out <- lr_test_pareto(x, xm = 1)
  expect_true(all(out$preferred == "pareto"))
})

test_that("lr_test_pareto works with a subset of alternatives", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  out <- lr_test_pareto(x, xm = 1, alternatives = c("exponential", "lognormal"))
  expect_equal(nrow(out), 2)
})

test_that("lr_test_pareto p_values are in [0, 1]", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  out <- lr_test_pareto(x, xm = 1)
  expect_true(all(out$p_value >= 0 & out$p_value <= 1))
})

test_that("lr_test_pareto validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(lr_test_pareto(as.character(x), xm = 1))
  expect_error(lr_test_pareto(c(1), xm = 1))
  expect_error(lr_test_pareto(x, xm = -1))
  expect_error(lr_test_pareto(na_data, xm = 1, na.rm = FALSE))
  expect_error(lr_test_pareto(x, xm = 1, alternatives = "poisson"))
})
