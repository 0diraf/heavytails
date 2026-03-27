##### Hill Estimator ####

#' Hill Estimator
#'
#' Hill estimator used to calculate the tail index (alpha) of input data.
#'
#' \deqn{\hat \alpha_H = \frac{1}{\frac{1}{k} \sum_{i=1}^{k} log(\frac{X_{(i)}}{X_{(k+1)}})}}
#'
#' where \eqn{X_{(1)} \ge X_{(2)} \ge \dots \ge X_{(n)}} are the order statistics
#' of the data (descending order).
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param k An integer specifying the number of top order statistics to use
#'   (the size of the tail). Must be strictly less than the sample size.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A single numeric scalar: Hill estimator calculation of the tail index \eqn{\alpha}.
#' @export
#'
#' @examples
#'
#' xmin <- 1
#' alpha <- 2
#' r <- runif(800, 0, 1)
#' x <- (xmin * r^(-1/(alpha)))
#' hill <- hill_estimator(data = x, k = 5)
#'
#'
#' @references
#'
#' Hill, B. M. (1975). A Simple General Approach to Inference About the Tail of a Distribution. \emph{The Annals of Statistics}, \bold{3}(5), 1163–1174. \url{http://www.jstor.org/stable/2958370}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 203-205) \doi{10.1017/9781009053730}
#'
hill_estimator <- function(data, k, na.rm = FALSE) {


  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(k) || length(k)!=1) {
    stop("`k` must be a numeric scalar")
  }

  if (length(data) <= 1 ) {
    stop("`data` must be a vector with length > 1")
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

  if (k < 1 || k >= length(data) - 1) {
    stop('k is out of valid range')
  }

  if (k==1) {
    warning('`hill_estimator` is undefined at k = 1')
  }

  n <- length(data)

  data <- sort(data, decreasing=T)

  Xi <- data[1:k]
  Xk1 <- data[k+1]

  hill <- 1/(mean(log(Xi)-log(Xk1)))

  return(hill)

}
