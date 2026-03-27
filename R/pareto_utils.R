#### Pareto Distribution Utilities ####

#' Generate Pareto Random Variates
#'
#' Generates \eqn{n} random samples from a Pareto(\eqn{x_m}, \eqn{\alpha})
#' distribution via inverse CDF: \eqn{x_m \cdot U^{-1/\alpha}} where
#' \eqn{U \sim \text{Uniform}(0, 1)}.
#'
#' @param n A non-negative integer: number of samples to generate.
#' @param alpha A positive numeric scalar: tail index.
#' @param xm A positive numeric scalar: scale parameter (lower bound).
#'
#' @returns A numeric vector of length \code{n}.
#' @export
#'
#' @examples
#'
#' x <- rpareto(n = 500, alpha = 2, xm = 1)
#'
#' @references
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy
#' Tails: Properties, Emergence, and Estimation}. Cambridge University Press.
#' \doi{10.1017/9781009053730}
#'
rpareto <- function(n, alpha, xm) {
  if (!is.numeric(n) || length(n) != 1 || n < 0 || n != floor(n)) {
    stop("`n` must be a non-negative integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a positive numeric scalar.")
  }
  if (!is.numeric(xm) || length(xm) != 1 || xm <= 0) {
    stop("`xm` must be a positive numeric scalar.")
  }
  if (n == 0) {
    return(numeric(0))
  }
  r <- stats::runif(n, 0, 1)
  return(xm * r^(-1 / alpha))
}

#' Pareto CDF
#'
#' Computes the cumulative distribution function of the Pareto(\eqn{x_m},
#' \eqn{\alpha}) distribution:
#' \deqn{F(x) = 1 - \left(\frac{x}{x_m}\right)^{-\alpha}, \quad x \ge x_m}
#'
#' @param x A numeric vector of quantiles.
#' @param xmin A positive numeric scalar: scale parameter (lower bound).
#' @param alpha A positive numeric scalar: tail index.
#'
#' @returns A numeric vector of CDF values in \eqn{[0, 1]}.
#' @export
#'
#' @examples
#'
#' pareto_cdf(x = c(1, 2, 5), xmin = 1, alpha = 2)
#'
pareto_cdf <- function(x, xmin, alpha) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }
  if (!is.numeric(xmin) || length(xmin) != 1 || xmin <= 0) {
    stop("`xmin` must be a positive numeric scalar.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a positive numeric scalar.")
  }
  ifelse(x < xmin, 0, 1 - (x / xmin)^(-alpha))
}

#' Pareto Density
#'
#' Computes the probability density function of the Pareto(\eqn{x_m},
#' \eqn{\alpha}) distribution:
#' \deqn{f(x) = \frac{\alpha \, x_m^{\alpha}}{x^{\alpha + 1}}, \quad x \ge x_m}
#'
#' @param x A numeric vector of quantiles.
#' @param alpha A positive numeric scalar: tail index.
#' @param xm A positive numeric scalar: scale parameter (lower bound).
#'
#' @returns A numeric vector of density values (zero for \eqn{x < x_m}).
#' @export
#'
#' @examples
#'
#' dpareto(x = c(1, 2, 5), alpha = 2, xm = 1)
#'
dpareto <- function(x, alpha, xm) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a positive numeric scalar.")
  }
  if (!is.numeric(xm) || length(xm) != 1 || xm <= 0) {
    stop("`xm` must be a positive numeric scalar.")
  }
  ifelse(x < xm, 0, alpha * xm^alpha / x^(alpha + 1))
}
