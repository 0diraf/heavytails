test_that("pot_estimator returns numeric scalar", {

  set.seed(1)
  x <- rweibull(n=1000, shape = .8, scale = 1)

  thr <- quantile(x, 0.95)

  est <- pot_estimator(x, u = thr)

  expect_type(est, "double")
  expect_length(est, 2)
})

test_that("pot_estimator errors for bad threshold", {

  x <- rweibull(n=500, shape = .8, scale = 1)

  expect_error(pot_estimator(x, u = -1))
  expect_error(pot_estimator(x, u = 10^6))  # above all data
})

test_that("pot_estimator behaves reasonably on Pareto data", {

  set.seed(1)
  x <- heavytails:::rpareto(1000, xm = 1, alpha = 2.0)

  true_xi <- 1/2

  thr <- quantile(x, 0.90)

  est <- pot_estimator(x, u = thr)[1]

  expect_true(abs(est - true_xi) < 0.2)
})

test_that("pot_estimator validates input arguments", {
  x <- c(1)
  y <- 1
  z <- c(1,2,NA,NA,9,14,20,78)
  z1 <- c(1,NA,NA,NA)
  a <- c(1,4,2,4,5,1,9,10,18,72,88,129)
  threshold <- max(a)
  xi <- c(2,4,5)
  beta <- c(30, 8, 9)

  expect_error(pot_estimator(x, u = 1)) #Length of data must be greater than 1
  expect_error(pot_estimator(y, u = 1)) #Data must be a vector
  expect_error(pot_estimator(z, u = 3, na.rm=FALSE))  #NA values not handled
  expect_error(pot_estimator(z1, u = 1, na.rm=TRUE))  #Removal of NA results in length(z1) == 1.
  expect_error(pot_estimator(a, u = threshold))       #Threshold cannot be >= max(data)
  expect_error(pot_estimator(a, u = 9, start_xi = xi, start_beta = beta))   #Beta and xi must be numeric scalers

})
