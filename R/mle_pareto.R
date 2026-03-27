#### Parametric MLE for the Pareto Distribution ####

#' Parametric MLE for the Pareto Distribution
#'
#' Estimates the tail index \eqn{\alpha} of a Pareto(\eqn{x_m}, \eqn{\alpha})
#' distribution via maximum likelihood (Theorem 8.1 of Nair et al.).
#'
#' The MLE is:
#' \deqn{\hat{\alpha} = \frac{n}{\sum_{i=1}^{n} \log(X_i / x_m)}}
#'
#' Unlike the Hill estimator (which uses only the top \eqn{k} order statistics),
#' this estimator uses all \eqn{n} observations and assumes the entire sample
#' follows a Pareto distribution with known lower bound \eqn{x_m}.
#'
#' A finite-sample bias-corrected version (§8.3) uses \eqn{n - 1} in the
#' numerator:
#' \deqn{\hat{\alpha}^* = \frac{n - 1}{\sum_{i=1}^{n} \log(X_i / x_m)}}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param xm Optional positive numeric scalar. Lower bound of the Pareto
#'   support. If \code{NULL} (default), \code{min(data)} is used.
#' @param bias_corrected Logical. If \code{TRUE} (default), applies the
#'   finite-sample bias correction described in §8.3.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A named list with elements:
#' \itemize{
#'   \item{\code{alpha}:} Estimated tail index.
#'   \item{\code{xm}:} The lower bound used.
#'   \item{\code{n}:} Number of observations used (those \eqn{\ge x_m}).
#'   \item{\code{bias_corrected}:} Logical indicating whether bias correction
#'     was applied.
#' }
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(n = 1000, alpha = 2, xm = 1)
#' mle_pareto(x)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' (pp. 162-167) \doi{10.1017/9781009053730}
#'
mle_pareto <- function(data, xm = NULL, bias_corrected = TRUE, na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(xm) && (!is.numeric(xm) || length(xm) != 1 || xm <= 0)) {
    stop("`xm` must be a positive numeric scalar or NULL.")
  }

  if (!is.logical(bias_corrected) || length(bias_corrected) != 1) {
    stop("`bias_corrected` must be a logical scalar.")
  }

  if (length(data) <= 1) {
    stop("`data` must be a vector with length > 1")
  }

  if (anyNA(data) && na.rm == TRUE) {
    data <- data[!is.na(data)]
    if (length(data) <= 1) {
      stop("Removing NAs resulted in a data vector with length <= 1. Data must be a vector with length > 1")
    }
  }

  if (anyNA(data)) {
    stop("`data` must not contain NA values. Please remove NAs, or set na.rm = TRUE in the function call")
  }

  if (is.null(xm)) {
    xm <- min(data)
  }

  x_valid <- data[data >= xm]

  if (length(x_valid) == 0) {
    stop("No data values >= `xm`. Check that `xm` is within the data range.")
  }

  if (any(data < xm)) {
    warning(paste0(sum(data < xm), " observation(s) below `xm` were excluded."))
  }

  n <- length(x_valid)

  log_ratios <- log(x_valid / xm)
  sum_log <- sum(log_ratios)

  if (sum_log <= 0) {
    stop("Sum of log(data/xm) is non-positive. Check that `xm` is a valid lower bound.")
  }

  if (bias_corrected) {
    alpha_hat <- (n - 1) / sum_log
  } else {
    alpha_hat <- n / sum_log
  }

  result <- list(alpha = alpha_hat, xm = xm, n = n, bias_corrected = bias_corrected)
  class(result) <- "heavytails_mle"
  result
}

#' @export
print.heavytails_mle <- function(x, ...) {
  cat("Pareto MLE\n")
  cat("  alpha:           ", round(x$alpha, 4), "\n")
  cat("  xm:              ", round(x$xm, 4), "\n")
  cat("  n (used):        ", x$n, "\n")
  cat("  bias_corrected:  ", x$bias_corrected, "\n")
  invisible(x)
}
