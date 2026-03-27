test_that("pickands_estimator returns numeric scalar", {

  set.seed(1)
  x <- heavytails:::rpareto(2000, alpha = 2.5, xm = 1)

  est <- pickands_estimator(x, k = 100)

  expect_type(est, "double")
  expect_length(est, 1)
})

test_that("pickands_estimator k must be valid", {

  x <- heavytails:::rpareto(2000, alpha = 2.5, xm = 1)
  expect_error(pickands_estimator(x, k = 0))   # too small
  expect_error(pickands_estimator(x, k = 6000)) # too large
})

test_that("pickands_estimator is roughly correct on Pareto data", {

  set.seed(2)
  x <- heavytails:::rpareto(2000, alpha = 2.5, xm = 1)

  est <- pickands_estimator(x, k = 200)

  true_xi <- 1 / 2.5   # 0.4

  expect_true(abs(est - true_xi) < 0.2)
})

test_that("pickands_estimator validates input arguments", {
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)
  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)

  expect_error(pickands_estimator(x, k = 1)) #Length of data must be greater than 4
  expect_error(pickands_estimator(y, k = 1)) #Data must be a vector
  expect_error(pickands_estimator(z, k = 3, na.rm=FALSE))  #NA values not handled
  expect_error(pickands_estimator(z1, k = 1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.
  expect_error(pickands_estimator(a, k = 0))               #k < 1 not allowed

})
