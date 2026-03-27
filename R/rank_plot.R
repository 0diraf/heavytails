#### Rank Plot (Log-Log) ####

#' Rank Plot (Log-Log CCDF)
#'
#' Plots the empirical complementary CDF (CCDF) of the data on a log-log scale.
#' A power-law distribution appears as a straight line on this plot. If a fitted
#' \code{plfit()} result is supplied, the theoretical Pareto CCDF is overlaid.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param fit Optional. A list returned by \code{\link{plfit}}, used to overlay
#'   the fitted Pareto line. Must contain \code{alpha_hat} and \code{xmin_hat}.
#' @param log_scale Logical. If \code{TRUE} (default), axes are log-transformed
#'   (i.e., \eqn{\log x} vs \eqn{\log \hat{F}^c}).
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns A \code{data.frame} with columns \code{x} and \code{ccdf},
#'   returned invisibly.
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(800, alpha = 2, xm = 1)
#' fit <- plfit(x)
#' rank_plot(x, fit = fit)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' (pp. 176-179) \doi{10.1017/9781009053730}
#'
rank_plot <- function(data, fit = NULL, log_scale = TRUE, na.rm = FALSE, ...) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(fit) && (!is.list(fit) || is.null(fit$alpha_hat) || is.null(fit$xmin_hat))) {
    stop("`fit` must be a list with elements `alpha_hat` and `xmin_hat` (as returned by plfit()).")
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

  data <- data[data > 0]
  if (length(data) == 0) {
    stop("`data` has no positive values.")
  }

  x <- sort(data, decreasing = TRUE)
  n <- length(x)
  ccdf <- seq_len(n) / n   # i/n — rank-based CCDF (Eq. 8.3)

  if (log_scale) {
    plot_x <- log(x)
    plot_y <- log(ccdf)
    xlab_default <- "log(x)"
    ylab_default <- expression(log(hat(F)^c))
  } else {
    plot_x <- x
    plot_y <- ccdf
    xlab_default <- "x"
    ylab_default <- expression(hat(F)^c)
  }

  args <- list(x = plot_x, y = plot_y, type = "p", pch = 20, cex = 0.4,
               xlab = xlab_default, ylab = ylab_default)
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(graphics::plot, args)

  if (!is.null(fit)) {
    alpha_hat <- fit$alpha_hat
    xmin_hat  <- fit$xmin_hat
    x_tail <- x[x >= xmin_hat]
    theoretical_ccdf <- (x_tail / xmin_hat)^(-alpha_hat)
    if (log_scale) {
      graphics::lines(log(x_tail), log(theoretical_ccdf), col = "red", lwd = 1.5)
    } else {
      graphics::lines(x_tail, theoretical_ccdf, col = "red", lwd = 1.5)
    }
  }

  result <- data.frame(x = x, ccdf = ccdf)
  invisible(result)
}
