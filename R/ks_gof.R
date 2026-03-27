#### Goodness-of-Fit via Parametric Bootstrap KS Test ####

#' Goodness-of-Fit Test for the Pareto Distribution
#'
#' Tests whether a Pareto(\eqn{x_m}, \eqn{\alpha}) distribution is a good fit
#' for the data by computing a bootstrap p-value for the Kolmogorov-Smirnov
#' (KS) statistic (Step 2 of the Clauset et al. pipeline, §8.5).
#'
#' The p-value is the fraction of bootstrap KS statistics that exceed the
#' observed KS statistic. A large p-value (e.g., > 0.1) means the Pareto
#' hypothesis cannot be rejected.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param alpha A positive numeric scalar: the Pareto tail index. Typically
#'   obtained from \code{\link{mle_pareto}} or \code{\link{plfit}}.
#' @param xm A positive numeric scalar: the lower bound. Only
#'   \code{data[data >= xm]} is used.
#' @param n_boot A positive integer: number of bootstrap replicates. Default
#'   \code{1000}.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A named list with elements:
#' \itemize{
#'   \item{\code{ks_statistic}:} Observed KS distance.
#'   \item{\code{p_value}:} Bootstrap p-value.
#'   \item{\code{n_boot}:} Number of bootstrap replicates used.
#'   \item{\code{n}:} Number of observations used (those \eqn{\ge x_m}).
#' }
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(n = 500, alpha = 2, xm = 1)
#' fit <- mle_pareto(x)
#' ks_gof(x, alpha = fit$alpha, xm = fit$xm, n_boot = 100)
#'
#' @references
#'
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions
#' in empirical data. \emph{SIAM Review}, \bold{51}(4), 661-703.
#' \doi{10.1137/070710111}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' (pp. 194-196) \doi{10.1017/9781009053730}
#'
ks_gof <- function(data, alpha, xm, n_boot = 1000, na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a positive numeric scalar.")
  }

  if (!is.numeric(xm) || length(xm) != 1 || xm <= 0) {
    stop("`xm` must be a positive numeric scalar.")
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) {
    stop("`n_boot` must be a positive integer.")
  }

  n_boot <- as.integer(n_boot)

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

  x_tail <- data[data >= xm]
  n <- length(x_tail)

  if (n < 2) {
    stop("Not enough observations >= `xm` to perform the test.")
  }

  observed_ks <- .compute_ks(x_tail, alpha, xm)

  boot_ks <- vapply(seq_len(n_boot), function(j) {
    x_boot <- rpareto(n, alpha, xm)
    fit_boot <- tryCatch(
      mle_pareto(x_boot, xm = xm, bias_corrected = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit_boot)) return(NA_real_)
    .compute_ks(x_boot, fit_boot$alpha, xm)
  }, numeric(1))

  boot_ks <- boot_ks[!is.na(boot_ks)]
  p_value <- mean(boot_ks >= observed_ks)

  list(ks_statistic = observed_ks, p_value = p_value, n_boot = n_boot, n = n)
}

# Internal helper: KS distance between data and Pareto CDF
#' @keywords internal
.compute_ks <- function(x, alpha, xm) {
  x_sorted <- sort(x[x >= xm])
  n <- length(x_sorted)
  if (n == 0) return(NA_real_)
  ecdf_vals <- seq_len(n) / n
  theoretical_vals <- pareto_cdf(x_sorted, xmin = xm, alpha = alpha)
  max(abs(ecdf_vals - theoretical_vals))
}
