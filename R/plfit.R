#### PLFIT ####

#' Power-law fit (PLFIT) Algorithm
#'
#' This function implements the PLFIT algorithm as described by *Clauset et al.* to determine the value of \eqn{\hat k}. It minimizes the Kolmorogorov-Smirnoff (KS) distance between the empirical cumulative distribution function and the fitted power law.
#'
#' \deqn{D_{n,k} := \sup_{y \ge 1} |\frac{1}{k-1} \sum_{i=1}^{k-1} I (\frac{X_{(i)}}{X_{(k)}} > y) - y^{-\hat{\alpha}_{n,k}^H}|}
#'
#'  The above equation, as described by *Nair et al.*, is implemented in this function with an Empirical CDF instead of the empirical survival function, which is mathematical equivalent since they are both complements of each other.
#'
#' \deqn{D_{n,k} :=
#'\sup_{y \ge 1}
#'|
#'  \underbrace{
#'    \frac{1}{k-1}
#'    \sum_{i=1}^{k-1}
#'    I(\frac{X_{(i)}}{X_{(k)}} \le y)
#'  }_{\text{Empirical CDF}}
#'-
#'  \underbrace{
#'    (1 - y^{-\hat{\alpha}_{n,k}})
#'  }_{\text{Theoretical CDF}}|}
#'
#'  \deqn{\hat k = \text{argmin} (D_{n,k})}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param kmax Maximum number of top-order statistics. If kmax = -1, then kmax=(n-1) where n is the length of dataset
#' @param kmin Minimum number of top-order statistics to start with
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A named list containing the results of the PLFIT algorithm:
#' \itemize{
#'   \item{\code{k_hat}:} The optimal number of top-order statistics \eqn{\hat{k}}.
#'   \item{\code{alpha_hat}:} The estimated power-law exponent \eqn{\hat{\alpha}} corresponding to \eqn{\hat{k}}.
#'   \item{\code{xmin_hat}:} The minimum value \eqn{x_{\min} = X_{(\hat{k})}} above which the power law is fitted.
#'   \item{\code{ks_distance}:} The minimum Kolmogorov-Smirnov distance \eqn{D_{n,k}} found.
#' }
#'
#' @export
#'
#' @examples
#'
#' xmin <- 1
#' alpha <- 2
#' r <- runif(800, 0, 1)
#' x <- (xmin * r^(-1/(alpha)))
#' plfit_values <- plfit(data = x, kmax = -1, kmin = 2)
#'
#' @references
#'
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions in empirical data. \emph{SIAM Review}, \bold{51}(4), 661-703. \doi{10.1137/070710111}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 227-229) \doi{10.1017/9781009053730}
#'
plfit <- function(data, kmax=-1, kmin=2, na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(kmin) || length(kmin) != 1 || kmin < 2) {
    stop("`kmin` must be a numeric scalar >= 2.")
  }

  if (!is.numeric(kmax) || length(kmax) != 1) {
    stop("`kmax` must be a numeric scalar.")
  }

  if (length(data) <= 1) {
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

  if (all(data <= 0)) {
    stop("`data` must contain positive values for power-law fitting.")
  }

  if (any(data <= 0)) {
    warning("`data` has negative values. They will be excluded from calculations.")
    data <- data[data > 0]
  }

  n <- length(data)

  if (n < 3) {
    stop("Not enough data to fit.")
  }

  if (kmax == -1) {
    kmax <- n - 1
  }

  result <- .ks_xmin_core(data, kmin = kmin, kmax = kmax)

  return(list(
    k_hat = result$k_hat,
    alpha_hat = result$alpha_hat,
    xmin_hat = result$xmin_hat,
    ks_distance = result$ks_distance
  ))
}

