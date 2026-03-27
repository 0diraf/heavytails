test_that("moments_estimator returns numeric scalar", {

  set.seed(1)
  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)

  est <- moments_estimator(x, k = 100)

  expect_type(est, "double")
  expect_length(est, 1)
})


test_that("moments_estimator works with pareto dist", {
  set.seed(1)

  # Pareto tail (xi = 1/alpha)
  x_pareto <- heavytails:::rpareto(1000, alpha = 2, xm = 1)  # xi = 0.5
  est_pareto <- moments_estimator(x_pareto, k = 200)
  expect_true(abs(est_pareto - 0.5) < 0.2)

})

test_that("moments_estimator returns NA for degenerate top k", {
  x_ident <- rep(10, 20)
  expect_true(is.na(moments_estimator(x_ident, k = 10)))
})

test_that("moments_estimator increases with heavier tails", {
  x_light <- heavytails:::rpareto(1000, alpha = 5, xm = 1)  # xi = 0.2
  x_heavy <- heavytails:::rpareto(1000, alpha = 2, xm = 1)  # xi = 0.5
  est_light <- moments_estimator(x_light, k = 200)
  est_heavy <- moments_estimator(x_heavy, k = 200)
  expect_true(est_heavy > est_light)
})

test_that("moments_estimator errors for invalid k", {

  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)
  expect_error(moments_estimator(x, k = 0))
  expect_error(moments_estimator(x, k = 3000))
})

test_that("moments_estimator validates input arguments", {
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)
  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)

  expect_error(moments_estimator(x, k = 1))  #Length of data must be greater than 1
  expect_error(moments_estimator(y, k = 1))  #Data must be a vector
  expect_error(moments_estimator(z, k = 3, na.rm=FALSE))  #NA values not handled
  expect_error(moments_estimator(z1, k = 1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.
  expect_error(moments_estimator(a, k = 0))               #k < 1 not allowed

})
