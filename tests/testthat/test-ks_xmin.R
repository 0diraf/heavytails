test_that("ks_xmin returns a list with the correct fields", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  out <- ks_xmin(x)

  expect_type(out, "list")
  expect_true(all(c("xm", "ks_distance", "k_hat") %in% names(out)))
  expect_true(is.numeric(out$xm))
  expect_true(is.numeric(out$ks_distance))
  expect_true(is.numeric(out$k_hat))
})

test_that("ks_xmin and plfit return consistent results", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  xmin_out <- ks_xmin(x)
  plfit_out <- plfit(x)

  expect_equal(xmin_out$xm,          plfit_out$xmin_hat)
  expect_equal(xmin_out$k_hat,        plfit_out$k_hat)
  expect_equal(xmin_out$ks_distance,  plfit_out$ks_distance)
})

test_that("ks_xmin xm estimate is near true xm on Pareto data", {
  set.seed(1)
  x <- rpareto(1000, alpha = 2, xm = 1)
  out <- ks_xmin(x)
  expect_true(out$xm >= 1)       # xm_hat must be >= true xm
  expect_true(out$xm < 1.5)      # and reasonably close
})

test_that("ks_xmin validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(ks_xmin(as.character(x)))
  expect_error(ks_xmin(c(1)))
  expect_error(ks_xmin(x, kmin = 1))
  expect_error(ks_xmin(na_data, na.rm = FALSE))
})
