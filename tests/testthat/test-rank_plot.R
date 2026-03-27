test_that("rank_plot returns a data.frame with x and ccdf columns", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- rank_plot(x)
  dev.off(); unlink(tmp)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("x", "ccdf") %in% names(out)))
  expect_equal(nrow(out), 300)
  expect_true(all(out$ccdf > 0 & out$ccdf <= 1))
})

test_that("rank_plot works with and without a plfit overlay", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  fit <- plfit(x)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_no_error(rank_plot(x))
  expect_no_error(rank_plot(x, fit = fit))
  dev.off(); unlink(tmp)
})

test_that("rank_plot works with log_scale = FALSE", {
  set.seed(1)
  x <- rpareto(200, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_no_error(rank_plot(x, log_scale = FALSE))
  dev.off(); unlink(tmp)
})

test_that("rank_plot validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(rank_plot(as.character(x)))
  expect_error(rank_plot(c(1)))
  expect_error(rank_plot(na_data, na.rm = FALSE))
  expect_error(rank_plot(x, fit = list(alpha_hat = 2)))  # missing xmin_hat
})
