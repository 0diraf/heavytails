test_that("doublebootstrap runs on simple data", {
  set.seed(1)
  x <- rexp(1000)

  out <- doublebootstrap(x, r = 10, k_max_prop = 0.5)

  expect_type(out, "list")
  expect_true("k" %in% names(out))
  expect_true("alpha" %in% names(out))
  expect_true(is.numeric(out$k))
  expect_true(is.numeric(out$alpha))
})

test_that("k_hat is always between 2 and n-2", {
  set.seed(1)
  x <- rexp(1000)

  out <- doublebootstrap(x, r = 10)
  expect_true(out$k >= 2)
  expect_true(out$k <= length(x) - 2)
})

test_that("alpha estimate is positive", {
  set.seed(1)
  x <- rexp(1000)

  out <- doublebootstrap(x, r = 10)
  expect_true(out$alpha > 0)
})

test_that("doublebootstrap validates input arguments", {

  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)
  b <- c(0, 0, 0, 0, 0, 0, 0)
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)

  expect_error(doublebootstrap(a, n1 = 6, n2 = 6))        #n2 must be <= n1
  expect_error(doublebootstrap(a, k_max_prop = -1))       #'k_max_prop' must be in (0,1].
  expect_error(doublebootstrap(x)) #Length of data must be greater than or equal to 3
  expect_error(doublebootstrap(y)) #Data must be a vector
  expect_error(doublebootstrap(z, na.rm=FALSE))  #NA values not handled
  expect_error(doublebootstrap(z1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.
})



