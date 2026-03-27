test_that("moments_plot returns a data.frame with k and xi_hat columns", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  out <- moments_plot(x)
  dev.off(); unlink(tmp)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("k", "xi_hat") %in% names(out)))
  expect_true(nrow(out) > 0)
})

test_that("moments_plot works with and without xi_true", {
  set.seed(1)
  x <- rpareto(300, alpha = 2, xm = 1)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_no_error(moments_plot(x))
  expect_no_error(moments_plot(x, xi_true = 0.5))
  dev.off(); unlink(tmp)
})

test_that("moments_plot validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(moments_plot(as.character(x)))
  expect_error(moments_plot(c(1)))
  expect_error(moments_plot(na_data, na.rm = FALSE))
  expect_error(moments_plot(x, xi_true = c(0.5, 1)))
})
