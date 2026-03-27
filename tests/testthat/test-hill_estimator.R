test_that("hill_estimator returns numeric scalar", {
  set.seed(1)
  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)
  est <- hill_estimator(x, k = 50)

  expect_type(est, "double")
  expect_length(est, 1)
})

test_that("hill_estimator errors for invalid k", {

  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)
  expect_error(hill_estimator(x, k = 0))
  expect_error(hill_estimator(x, k = 3000))  # k too large
})


test_that("hill_estimator is stable for large samples", {
  set.seed(1)
  x <- heavytails:::rpareto(n = 1000, alpha = 2, xm = 1)

  est1 <- hill_estimator(x, k = 200)
  est2 <- hill_estimator(x, k = 300)

  expect_true(abs(est1 - est2) < 0.2)
})

test_that("hill_estimator validates input arguments", {
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)
  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)

  expect_error(hill_estimator(x, k = 1)) #Length of data must be greater than 1
  expect_error(hill_estimator(y, k = 1)) #Data must be a vector
  expect_error(hill_estimator(z, k = 3, na.rm=FALSE))  #NA values not handled
  expect_error(hill_estimator(z1, k = 1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.
  expect_error(hill_estimator(a, k = 0))               #k < 1 not allowed

})
