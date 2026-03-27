#### Pareto QQ Plot ####

#' Pareto QQ Plot
#'
#' Produces a QQ plot comparing the empirical quantiles of the data (filtered
#' to \eqn{x \ge x_m}) against the theoretical quantiles of a
#' Pareto(\eqn{x_m}, \eqn{\alpha}) distribution. Points falling close to the
#' 45-degree reference line indicate a good Pareto fit.
#'
#' The theoretical quantile for the \eqn{i}-th order statistic is:
#' \deqn{q_i = x_m \left(\frac{n - i + 1}{n + 1}\right)^{-1/\alpha}}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param alpha A positive numeric scalar: the Pareto tail index (as returned
#'   by \code{\link{hill_estimator}} or \code{\link{mle_pareto}}).
#' @param xm Optional numeric scalar. Lower threshold; only data \eqn{\ge x_m}
#'   are used. If \code{NULL} (default), \code{min(data)} is used.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns A \code{data.frame} with columns \code{empirical} and
#'   \code{theoretical}, returned invisibly.
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(800, alpha = 2, xm = 1)
#' qq_pareto(x, alpha = 2, xm = 1)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' (pp. 191-194) \doi{10.1017/9781009053730}
#'
qq_pareto <- function(data, alpha, xm = NULL, na.rm = FALSE, ...) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a positive numeric scalar.")
  }

  if (!is.null(xm) && (!is.numeric(xm) || length(xm) != 1 || xm <= 0)) {
    stop("`xm` must be a positive numeric scalar or NULL.")
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

  x_filtered <- data[data >= xm]

  if (length(x_filtered) == 0) {
    stop("No data values >= `xm`.")
  }

  x_sorted <- sort(x_filtered)
  n <- length(x_sorted)

  # Theoretical quantiles: inverse Pareto CDF (Eq. 8.10)
  i <- seq_len(n)
  q_theoretical <- xm * ((n - i + 1) / (n + 1))^(-1 / alpha)

  args <- list(x = q_theoretical, y = x_sorted, type = "p", pch = 20, cex = 0.5,
               xlab = "Theoretical quantiles", ylab = "Empirical quantiles")
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(graphics::plot, args)

  graphics::abline(0, 1, lty = 2)

  result <- data.frame(empirical = x_sorted, theoretical = q_theoretical)
  invisible(result)
}
