test_that("hill_plot returns a data.frame with k and alpha_hat columns", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- hill_plot(x)
  dev.off(); unlink(tmp)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("k", "alpha_hat") %in% names(out)))
  expect_true(nrow(out) > 0)
  expect_true(all(out$alpha_hat > 0))
})

test_that("hill_plot respects k_range argument", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  k_range <- 10:50
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- hill_plot(x, k_range = k_range)
  dev.off(); unlink(tmp)

  expect_equal(out$k, k_range)
})

test_that("hill_plot works with and without alpha_true", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_no_error(hill_plot(x))
  expect_no_error(hill_plot(x, alpha_true = 2))
  dev.off(); unlink(tmp)
})

test_that("hill_plot validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(hill_plot(as.character(x)))
  expect_error(hill_plot(c(1)))
  expect_error(hill_plot(na_data, na.rm = FALSE))
  expect_error(hill_plot(x, alpha_true = c(1, 2)))
})
