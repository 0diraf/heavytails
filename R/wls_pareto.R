#### Weighted Least Squares Pareto Estimator ####

#' Weighted Least Squares Estimator for the Pareto Tail Index
#'
#' Estimates the Pareto tail index \eqn{\alpha} via weighted least squares (WLS)
#' regression on the log-log rank plot (Theorem 8.5 of Nair et al.). The WLS
#' weights \eqn{w_i = 1 / \log(X_{(i)} / x_m)} downweight noisy tail
#' observations relative to OLS, recovering the MLE asymptotically.
#'
#' The WLS estimate is:
#' \deqn{\hat{\alpha}_{WLS} = -\frac{\sum_i w_i \log(\hat{F}^c_i) \log(X_{(i)}/x_m)}{\sum_i w_i (\log(X_{(i)}/x_m))^2}}
#'
#' If \code{plot = TRUE}, the rank plot is drawn with both the WLS and OLS
#' fitted lines, reproducing Figure 8.9 of Nair et al.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param xm Optional positive numeric scalar. Lower bound. If \code{NULL}
#'   (default), \code{min(data)} is used.
#' @param plot Logical. If \code{TRUE} (default), draws the log-log rank plot
#'   with fitted WLS and OLS lines.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param ... Additional graphical arguments passed to
#'   \code{\link[graphics]{plot}} (only used when \code{plot = TRUE}).
#'
#' @returns A named list with elements:
#' \itemize{
#'   \item{\code{alpha_wls}:} WLS estimate of the tail index.
#'   \item{\code{alpha_ols}:} OLS estimate (unweighted) for comparison.
#'   \item{\code{xm}:} The lower bound used.
#' }
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(n = 500, alpha = 2, xm = 1)
#' wls_pareto(x)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' (pp. 167-173) \doi{10.1017/9781009053730}
#'
wls_pareto <- function(data, xm = NULL, plot = TRUE, na.rm = FALSE, ...) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(xm) && (!is.numeric(xm) || length(xm) != 1 || xm <= 0)) {
    stop("`xm` must be a positive numeric scalar or NULL.")
  }

  if (!is.logical(plot) || length(plot) != 1) {
    stop("`plot` must be a logical scalar.")
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

  x_sorted <- sort(data[data >= xm])
  n <- length(x_sorted)

  if (n < 2) {
    stop("Not enough data values >= `xm` to compute WLS estimate.")
  }

  i <- seq_len(n)
  ccdf_i <- (n + 1 - i) / n
  log_x  <- log(x_sorted / xm)
  log_ccdf <- log(ccdf_i)

  # Points where log_x == 0 (x == xm) carry no slope information; exclude them
  valid <- log_x > 0
  lx <- log_x[valid]
  lc <- log_ccdf[valid]

  if (length(lx) < 2) {
    stop("Too few distinct values > `xm` for WLS estimation.")
  }

  w <- 1 / lx   # weights: w_i = 1 / log(x_i / xm)

  # Closed-form WLS slope (intercept forced through origin on log-log scale)
  alpha_wls <- -sum(w * lc * lx) / sum(w * lx^2)

  # OLS slope for comparison
  alpha_ols <- -sum(lc * lx) / sum(lx^2)

  if (plot) {
    args <- list(x = log_x, y = log_ccdf, type = "p", pch = 20, cex = 0.4,
                 xlab = expression(log(x / x[m])), ylab = expression(log(hat(F)^c)),
                 main = "Log-log rank plot: WLS vs OLS")
    extra <- list(...)
    args[names(extra)] <- extra
    do.call(graphics::plot, args)

    x_line <- c(0, max(log_x))
    graphics::lines(x_line, -alpha_wls * x_line, col = "red",  lwd = 1.5, lty = 1)
    graphics::lines(x_line, -alpha_ols  * x_line, col = "blue", lwd = 1.5, lty = 2)
    graphics::legend("topright",
                     legend = c(paste0("WLS: alpha=", round(alpha_wls, 3)),
                                paste0("OLS: alpha=", round(alpha_ols, 3))),
                     col = c("red", "blue"), lty = c(1, 2), lwd = 1.5,
                     bty = "n")
  }

  result <- list(alpha_wls = alpha_wls, alpha_ols = alpha_ols, xm = xm)
  class(result) <- "heavytails_wls"
  result
}

#' @export
print.heavytails_wls <- function(x, ...) {
  cat("Pareto WLS\n")
  cat("  alpha_wls:  ", round(x$alpha_wls, 4), "\n")
  cat("  alpha_ols:  ", round(x$alpha_ols, 4), "\n")
  cat("  xm:         ", round(x$xm, 4), "\n")
  invisible(x)
}
