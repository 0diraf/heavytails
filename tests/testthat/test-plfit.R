test_that("plfit runs and returns expected fields", {

  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)

  out <- plfit(x)

  expect_type(out, "list")
  expect_true("k_hat" %in% names(out))
  expect_true("alpha_hat"  %in% names(out))
  expect_true("xmin_hat" %in% names(out))
  expect_true("ks_distance"  %in% names(out))
  expect_true(is.numeric(out$k_hat))
  expect_true(is.numeric(out$alpha_hat))
  expect_true(is.numeric(out$xmin_hat))
  expect_true(is.numeric(out$ks_distance))
})

test_that("plfit determines alpha of synthetic data", {
  set.seed(1)
  x <- heavytails:::rpareto(1000, alpha = 2.5, xm = 1)

  out <- plfit(x)

  expect_true(abs(out$alpha_hat - 2.5) < 0.3)
})

test_that("plfit validates input arguments", {

  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)
  b <- c(0, 0, 0, 0, 0, 0, 0)
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)

  expect_error(plfit(a, kmin = 1))                        #kmin must be greater than or equal to 2
  expect_error(plfit(a, kmin = c(2,3), kmax = c(2,3,4)))  #kmin, kmax must be numeric scalars
  expect_error(plfit(b)) #data must not be all 0
  expect_error(plfit(x)) #Length of data must be greater than or equal to 3
  expect_error(plfit(y)) #Data must be a vector
  expect_error(plfit(z, na.rm=FALSE))  #NA values not handled
  expect_error(plfit(z1, k = 1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.

})

