test_that("mle_pareto returns a list with the correct fields", {
  set.seed(1)
  x <- rpareto(500, alpha = 2, xm = 1)
  out <- mle_pareto(x)

  expect_type(out, "list")
  expect_true(all(c("alpha", "xm", "n", "bias_corrected") %in% names(out)))
  expect_type(out$alpha, "double")
  expect_length(out$alpha, 1)
  expect_true(out$alpha > 0)
  expect_true(out$n == 500)
})

test_that("mle_pareto estimate is close to true alpha on large sample", {
  set.seed(1)
  x <- rpareto(5000, alpha = 2, xm = 1)
  out <- mle_pareto(x)
  expect_true(abs(out$alpha - 2) < 0.1)
})

test_that("mle_pareto bias_corrected=TRUE gives a different result than FALSE on small samples", {
  set.seed(1)
  x <- rpareto(20, alpha = 2, xm = 1)
  out_corrected   <- mle_pareto(x, bias_corrected = TRUE)
  out_uncorrected <- mle_pareto(x, bias_corrected = FALSE)
  expect_false(isTRUE(all.equal(out_corrected$alpha, out_uncorrected$alpha)))
})

test_that("mle_pareto respects supplied xm and warns on excluded values", {
  set.seed(1)
  x <- rpareto(200, alpha = 2, xm = 1)
  x_with_low <- c(x, 0.1, 0.5)
  expect_warning(mle_pareto(x_with_low, xm = 1))
})

test_that("mle_pareto has correct S3 class", {
  x <- rpareto(100, alpha = 2, xm = 1)
  out <- mle_pareto(x)
  expect_s3_class(out, "heavytails_mle")
})

test_that("mle_pareto validates inputs", {
  x <- rpareto(100, alpha = 2, xm = 1)
  na_data <- c(x[1:5], NA, x[6:10])

  expect_error(mle_pareto(as.character(x)))
  expect_error(mle_pareto(c(1)))
  expect_error(mle_pareto(x, xm = -1))
  expect_error(mle_pareto(na_data, na.rm = FALSE))
  expect_error(mle_pareto(x[1:3], xm = 999))  # all data < xm
})
