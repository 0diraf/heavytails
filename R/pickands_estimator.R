#### Pickands Estimator ####


#' Pickands Estimator
#'
#' Pickands estimator to calculate \eqn{\xi} for the input data.
#'
#'
#' \deqn{\hat{\xi}_{P}  = \frac{1}{\log 2} \log ( \frac{X_{(k)} - X_{(2k)}}{X_{(2k)} - X_{(4k)}})}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param k An integer specifying the number of top order statistics to use
#'   (the size of the tail). Must be strictly less than the sample size.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#'
#' @returns A single numeric scalar: Pickands estimator calculation of the shape parameter \eqn{\xi}.
#' @export
#'
#' @examples
#'
#' xmin <- 1
#' alpha <- 2
#' r <- runif(800, 0, 1)
#' x <- (xmin * r^(-1/(alpha)))
#' pickands <- pickands_estimator(data = x, k = 5)
#'
#' @references
#'
#' Pickands, J. (1975). Statistical Inference Using Extreme Order Statistics. \emph{The Annals of Statistics}, \bold{3}(1), 119–131. \url{http://www.jstor.org/stable/2958083}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 219-221) \doi{10.1017/9781009053730}
#'
pickands_estimator <- function(data, k, na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(k) || length(k)!=1) {
    stop("`k` must be a numeric scalar")
  }

  if (length(data) <= 4 ) {
    stop("`data` must be a vector with length > 4")
  }

  if (anyNA(data) && na.rm==TRUE) {
    data <- data[!is.na(data)]
    if (length(data) <= 1) {
      stop("Removing NAs resulted in a data vector with length <= 1. Data must be a vector with length > 1")
    }
  }

  if (anyNA(data)) {
    stop("`data` must not contain NA values. Please remove NAs, or set na.rm = TRUE in the function call")
  }

  n <- length(data)

  if (k<1 || 4*k>=n) {
    stop("k is out of valid range")
  }

  data <- sort(data, decreasing = T)
  Xk <- data[k]
  X2k <- data[2*k]
  X4k <- data[4*k]

  denominator <- X2k-X4k
  numerator <- Xk - X2k

  if (denominator < 0 || numerator < 0) {
    stop("Pickands estimator is undefined for non-positive tail spacing.")
  }

  if (denominator == 0) {
    stop("Division by 0 in the calculation of the Pickands estimator")
  }

  division <- numerator/denominator

  if (division <= 0) {
    stop("Log undefined: Ratio term became 0")
  }

  gamma_hat <- 1/log(2)*log(division)

  return(gamma_hat)

}

