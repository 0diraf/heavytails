#### KS-Based x_min Estimation ####

#' Estimate the Power-Law Lower Bound via KS Minimization
#'
#' Estimates the lower bound \eqn{\hat{x}_m} of a power-law regime by finding
#' the order statistic that minimizes the Kolmogorov-Smirnov distance between
#' the empirical distribution and the fitted Pareto (Step 1 of the Clauset
#' et al. pipeline, §8.5).
#'
#' This function extracts and exposes the core loop from \code{\link{plfit}},
#' allowing \eqn{\hat{x}_m} estimation as a standalone step — useful as input
#' to \code{\link{mle_pareto}}, \code{\link{wls_pareto}}, or
#' \code{\link{ks_gof}}.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param kmax Maximum number of top order statistics to consider. If
#'   \code{-1} (default), uses \code{n - 1}.
#' @param kmin Minimum number of top order statistics. Default \code{2}.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A named list with elements:
#' \itemize{
#'   \item{\code{xm}:} Estimated lower bound \eqn{\hat{x}_m = X_{(\hat{k})}}.
#'   \item{\code{ks_distance}:} Minimum KS distance achieved.
#'   \item{\code{k_hat}:} The optimal \eqn{\hat{k}}.
#' }
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(n = 500, alpha = 2, xm = 1)
#' ks_xmin(x)
#'
#' @references
#'
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions
#' in empirical data. \emph{SIAM Review}, \bold{51}(4), 661-703.
#' \doi{10.1137/070710111}
#'
ks_xmin <- function(data, kmax = -1, kmin = 2, na.rm = FALSE) {

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

  if (anyNA(data) && na.rm == TRUE) {
    data <- data[!is.na(data)]
    if (length(data) <= 1) {
      stop("Removing NAs resulted in a data vector with length <= 1. Data must be a vector with length > 1")
    }
  }

  if (anyNA(data)) {
    stop("`data` must not contain NA values. Please remove NAs, or set na.rm = TRUE in the function call")
  }

  if (all(data <= 0)) {
    stop("`data` must contain positive values.")
  }

  if (any(data <= 0)) {
    warning("`data` has non-positive values. They will be excluded from calculations.")
    data <- data[data > 0]
  }

  n <- length(data)

  if (n < 3) {
    stop("Not enough data to estimate x_min.")
  }

  if (kmax == -1) {
    kmax <- n - 1
  }

  result <- .ks_xmin_core(data, kmin = kmin, kmax = kmax)

  list(xm = result$xmin_hat, ks_distance = result$ks_distance, k_hat = result$k_hat)
}

# Internal core loop shared by ks_xmin() and plfit()
#' @keywords internal
.ks_xmin_core <- function(data, kmin, kmax) {
  x <- sort(data, decreasing = TRUE)
  n <- length(x)

  ks_distances <- rep(NA_real_, kmax)
  alphas        <- rep(NA_real_, kmax)

  for (k in kmin:kmax) {
    k1 <- k - 1
    current_xmin <- x[k]
    scaled <- x[1:k1] / current_xmin
    log_s  <- log(scaled)
    xi_est <- mean(log_s[is.finite(log_s)])
    current_alpha <- if (xi_est == 0) Inf else 1 / xi_est
    alphas[k] <- current_alpha

    ecdf_vals <- (1:k1) / k1
    scaled_sorted <- sort(scaled, decreasing = FALSE)
    theoretical_cdf_vals <- pareto_cdf(scaled_sorted, xmin = 1, alpha = current_alpha)
    ks_distances[k] <- max(abs(ecdf_vals - theoretical_cdf_vals))
  }

  k_hat <- which.min(ks_distances)

  list(k_hat = k_hat, alpha_hat = alphas[k_hat],
       xmin_hat = x[k_hat], ks_distance = ks_distances[k_hat])
}
