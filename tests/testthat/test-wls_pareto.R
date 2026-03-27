test_that("wls_pareto returns a list with the correct fields", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  out <- wls_pareto(x, plot = FALSE)

  expect_type(out, "list")
  expect_true(all(c("alpha_wls", "alpha_ols", "xm") %in% names(out)))
  expect_true(out$alpha_wls > 0)
  expect_true(out$alpha_ols > 0)
})

test_that("wls_pareto alpha_wls is close to mle_pareto on the same data", {
  set.seed(1)
  x <- rpareto(2000, alpha = 2, xm = 1)
  wls <- wls_pareto(x, plot = FALSE)
  mle <- mle_pareto(x, bias_corrected = FALSE)
  expect_true(abs(wls$alpha_wls - mle$alpha) < 0.05)
})

test_that("wls_pareto has correct S3 class", {
  x <- rpareto(100, alpha = 2, xm = 1)
  out <- wls_pareto(x, plot = FALSE)
  expect_s3_class(out, "heavytails_wls")
})

test_that("wls_pareto plot=TRUE runs without error", {
  set.seed(1)
  x <- rpareto(200, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_no_error(wls_pareto(x, plot = TRUE))
  dev.off()
  unlink(tmp)
})

test_that("wls_pareto validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(wls_pareto(as.character(x)))
  expect_error(wls_pareto(c(1)))
  expect_error(wls_pareto(x, xm = -1))
  expect_error(wls_pareto(na_data, na.rm = FALSE))
})
