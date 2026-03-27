#' Moments Estimator
#'
#' Moments estimator to calculate \eqn{\xi} for the input data.
#'
#'
#' \deqn{\hat \xi_{ME} = \underbrace{\hat \xi_{k,n}^{H,1}}_{T_1} + \underbrace{1 - \frac{1}{2}(1-\frac{(\hat \xi_{k,n}^{H,1})^2}{\hat \xi_{k,n}^{H,2}})^{-1}}_{T_2}}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param k An integer specifying the number of top order statistics to use
#'   (the size of the tail). Must be strictly less than the sample size.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param eps numeric, factor added to the denominator to avoid division by zero. Default value is 1e-12.
#'
#' @returns A single numeric scalar: Moments estimator calculation of the shape parameter \eqn{\xi}.
#' @export
#'
#' @examples
#'
#' xmin <- 1
#' alpha <- 2
#' r <- runif(800, 0, 1)
#' x <- (xmin * r^(-1/(alpha)))
#' moments <- moments_estimator(data = x, k = 5)
#'
#' @references
#'
#' Dekkers, A. L. M., Einmahl, J. H. J., & De Haan, L. (1989). A Moment Estimator for the Index of an Extreme-Value Distribution. \emph{The Annals of Statistics}, \bold{17}(4), 1833–1855. \url{http://www.jstor.org/stable/2241667}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 216-219) \doi{10.1017/9781009053730}
#'
moments_estimator <- function(data, k, na.rm = FALSE, eps = 1e-12) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(k) || length(k)!=1) {
    stop("`k` must be a numeric scalar")
  }

  if (!is.numeric(eps) || length(eps)!=1) {
    stop("`eps` must be a numeric scalar")
  }

  if (length(data) <= 1 ) {
    stop("Data must be a vector with length > 1")
  }

  if (anyNA(data) && na.rm==TRUE) {
    data <- data[!is.na(data)]
    if (length(data) <= 1) {
      stop("Removing NAs resulted in a data vector with length <= 1. Data must be a vector with length > 1")
    }
  }

  if (anyNA(data)) {
    stop("Data must not contain NA values. Please remove NAs, or set na.rm = TRUE in the function call")
  }


  n <- length(data)

  if (k < 1 || k >= n - 1) {
    stop('k is out of valid range')
  }

  data <- sort(data, decreasing=T)

  Xi <- data[1:k]
  Xk1 <- data[k+1]

  log_diff <- log(Xi) - log(Xk1)

  if (all(log_diff < eps)) {
    return(NA)
  }

  T1 <- mean(log_diff)
  squared <- mean(log_diff^2)

  ratio <- T1^(2)/(squared+eps)

  moments <- T1 + 1 - 0.5*(1-ratio)^(-1)

  return(moments)
}
