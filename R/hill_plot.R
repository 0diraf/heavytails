#### Hill Plot ####

#' Hill Plot
#'
#' Plots the Hill estimator of the tail index \eqn{\hat{\alpha}} as a function
#' of the number of top order statistics \eqn{k}. A stable plateau in this plot
#' is used to visually select a suitable value of \eqn{k}.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param k_range An integer vector specifying which values of \eqn{k} to
#'   evaluate. If \code{NULL} (default), uses \code{2:floor(length(data)/2)}.
#' @param alpha_true Optional numeric scalar. If supplied, a horizontal
#'   reference line at the true \eqn{\alpha} is added to the plot.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns A \code{data.frame} with columns \code{k} and \code{alpha_hat},
#'   returned invisibly. Users who prefer ggplot2 can capture this output
#'   and re-plot.
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(800, alpha = 2, xm = 1)
#' result <- hill_plot(x)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' \doi{10.1017/9781009053730}
#'
hill_plot <- function(data, k_range = NULL, alpha_true = NULL, na.rm = FALSE, ...) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(alpha_true) && (!is.numeric(alpha_true) || length(alpha_true) != 1)) {
    stop("`alpha_true` must be a numeric scalar or NULL.")
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

  n <- length(data)

  if (is.null(k_range)) {
    k_range <- 2:floor(n / 2)
  }

  if (!is.numeric(k_range) || length(k_range) < 1) {
    stop("`k_range` must be a non-empty numeric vector.")
  }

  # Clamp k_range so no value exceeds n - 2 (hill_estimator requires k < n - 1)
  k_range <- k_range[k_range >= 2 & k_range <= n - 2]

  if (length(k_range) == 0) {
    stop("`k_range` contains no valid k values after clamping to [2, n-2].")
  }

  alpha_hat <- sapply(k_range, function(k) hill_estimator(data, k))

  args <- list(x = k_range, y = alpha_hat, type = "l",
               xlab = "k", ylab = expression(hat(alpha)))
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(graphics::plot, args)

  if (!is.null(alpha_true)) {
    graphics::abline(h = alpha_true, lty = 2, col = "red")
  }

  result <- data.frame(k = k_range, alpha_hat = alpha_hat)
  invisible(result)
}
