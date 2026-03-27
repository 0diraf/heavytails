test_that("qq_pareto returns a data.frame with empirical and theoretical columns", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- qq_pareto(x, alpha = 2, xm = 1)
  dev.off(); unlink(tmp)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("empirical", "theoretical") %in% names(out)))
  expect_equal(nrow(out), 300)
  expect_true(all(out$empirical > 0))
  expect_true(all(out$theoretical > 0))
})

test_that("qq_pareto empirical and theoretical quantiles are close on Pareto data", {
  set.seed(1)
  x <- rpareto(1000, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- qq_pareto(x, alpha = 2, xm = 1)
  dev.off(); unlink(tmp)

  # On the log scale, the ratio should be close to 1 for most quantiles
  log_ratio <- log(out$empirical) - log(out$theoretical)
  expect_true(median(abs(log_ratio)) < 0.15)
})

test_that("qq_pareto uses min(data) as xm when xm is NULL", {
  set.seed(1)
  x <- rpareto(200, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out_null <- qq_pareto(x, alpha = 2, xm = NULL)
  out_min  <- qq_pareto(x, alpha = 2, xm = min(x))
  dev.off(); unlink(tmp)

  expect_equal(out_null, out_min)
})

test_that("qq_pareto validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(qq_pareto(as.character(x), alpha = 2))
  expect_error(qq_pareto(c(1), alpha = 2))
  expect_error(qq_pareto(x, alpha = -1))
  expect_error(qq_pareto(x, alpha = 2, xm = -1))
  expect_error(qq_pareto(na_data, alpha = 2, na.rm = FALSE))
})
