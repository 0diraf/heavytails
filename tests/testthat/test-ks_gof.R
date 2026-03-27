test_that("ks_gof returns a list with the correct fields", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  out <- ks_gof(x, alpha = 2, xm = 1, n_boot = 99)

  expect_type(out, "list")
  expect_true(all(c("ks_statistic", "p_value", "n_boot", "n") %in% names(out)))
  expect_true(is.numeric(out$ks_statistic))
  expect_true(out$p_value >= 0 && out$p_value <= 1)
  expect_equal(out$n_boot, 99L)
  expect_equal(out$n, 300L)
})

test_that("ks_gof gives a high p-value on genuine Pareto data", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  out <- ks_gof(x, alpha = 2, xm = 1, n_boot = 199)
  expect_true(out$p_value > 0.05)
})

test_that("ks_gof validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(ks_gof(as.character(x), alpha = 2, xm = 1))
  expect_error(ks_gof(c(1),            alpha = 2, xm = 1))
  expect_error(ks_gof(x, alpha = -1,   xm = 1))
  expect_error(ks_gof(x, alpha = 2,    xm = -1))
  expect_error(ks_gof(x, alpha = 2,    xm = 1, n_boot = 0))
  expect_error(ks_gof(na_data, alpha = 2, xm = 1, na.rm = FALSE))
})
