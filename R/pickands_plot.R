#### Pickands Plot ####

#' Pickands Estimator Plot
#'
#' Plots the Pickands estimator of the shape parameter \eqn{\hat{\xi}} as a
#' function of the number of top order statistics \eqn{k}. A stable plateau
#' indicates a suitable choice of \eqn{k}.
#'
#' The Pickands estimator requires \eqn{4k < n}, so the default \code{k_range}
#' upper bound is \code{floor(n/4) - 1}.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param k_range An integer vector specifying which values of \eqn{k} to
#'   evaluate. If \code{NULL} (default), uses \code{2:floor(length(data)/4 - 1)}.
#' @param xi_true Optional numeric scalar. If supplied, a horizontal reference
#'   line at the true \eqn{\xi} is added.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns A \code{data.frame} with columns \code{k} and \code{xi_hat},
#'   returned invisibly.
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(800, alpha = 2, xm = 1)
#' pickands_plot(x)
#'
#' @references
#'
#' Pickands, J. (1975). Statistical Inference Using Extreme Order Statistics.
#' \emph{The Annals of Statistics}, \bold{3}(1), 119–131.
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' \doi{10.1017/9781009053730}
#'
pickands_plot <- function(data, k_range = NULL, xi_true = NULL, na.rm = FALSE, ...) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(xi_true) && (!is.numeric(xi_true) || length(xi_true) != 1)) {
    stop("`xi_true` must be a numeric scalar or NULL.")
  }

  if (length(data) <= 4) {
    stop("`data` must be a vector with length > 4")
  }

  if (anyNA(data) && na.rm == TRUE) {
    data <- data[!is.na(data)]
    if (length(data) <= 4) {
      stop("Removing NAs resulted in a data vector with length <= 4. Data must be a vector with length > 4")
    }
  }

  if (anyNA(data)) {
    stop("`data` must not contain NA values. Please remove NAs, or set na.rm = TRUE in the function call")
  }

  n <- length(data)

  if (is.null(k_range)) {
    k_max <- floor(n / 4) - 1
    if (k_max < 2) stop("`data` is too short for the Pickands estimator (need n >= 12).")
    k_range <- 2:k_max
  }

  if (!is.numeric(k_range) || length(k_range) < 1) {
    stop("`k_range` must be a non-empty numeric vector.")
  }

  # Pickands requires 4*k < n
  k_range <- k_range[k_range >= 1 & 4 * k_range < n]

  if (length(k_range) == 0) {
    stop("`k_range` contains no valid k values satisfying 4*k < n.")
  }

  xi_hat <- sapply(k_range, function(k) {
    tryCatch(pickands_estimator(data, k), error = function(e) NA_real_)
  })

  args <- list(x = k_range, y = xi_hat, type = "l",
               xlab = "k", ylab = expression(hat(xi)))
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(graphics::plot, args)

  if (!is.null(xi_true)) {
    graphics::abline(h = xi_true, lty = 2, col = "red")
  }

  result <- data.frame(k = k_range, xi_hat = xi_hat)
  invisible(result)
}
