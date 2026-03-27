
#### Peaks over Threshold (POT) Estimator ####


#' Negative Log likelihood of Generalized Pareto Distribution
#'
#' Helper function for pot_estimator(). Returns the \eqn{\xi} and \eqn{\beta} that minimize the negative log-likelihood of the Generalized Pareto Distribution (GPD).
#'
#' \deqn{l(\xi,\beta)=-n\log(\beta) - (\frac{1}{\xi} + 1)\sum_{i=1}^{n}  \log(1 + \xi \frac{x_i}{\beta})}
#'
#' @param params Vector containing initial values of \eqn{\xi} and \eqn{\beta}
#' @param data Original dataset
#'
#' @returns Negative log-likelihood of the GPD.
#' @keywords internal
gpd_lg_likelihood <- function(params, data) {

  xi <- params[1]
  beta <- params[2]

  # NA, negative checks
  if (is.na(xi) || is.na(beta) || beta <=0 || (xi!=0 && any(1 + xi*data/beta<=0,
                                                            na.rm = TRUE))) {
    return(1e10)
  }

  n <-length(data)


  if (abs(xi) < 1e-6) {
    # Exponential
    log_lik <- -n*log(beta)-sum(data/beta)
  } else {
    # GPD
    log_lik <- -n*log(beta)-(1/xi+1)*sum(log(1+xi*data/beta))
  }

  return(-log_lik)
}

#' Peaks-over-threshold (POT) Estimator
#'
#'
#' This function chooses the \eqn{\hat{\xi}_{k}} and \eqn{\hat \beta} that minimize the negative log likelihood of the Generalized Pareto Distribution (GPD).
#'
#'The PDF of a excess data point \eqn{x_i} is given by:
#'
#' \deqn{f(x_i;\xi, \beta) = \frac{1}{\beta} \left(1 + \xi \frac{x_i}{\beta}\right)^{-\left(\frac{1}{\xi} + 1\right)}}
#'
#' If we apply \eqn{log} to the above equation we get:
#'
#' \deqn{l(x_i;\xi, \beta)=-\log(\beta) - (\frac{1}{\xi} + 1) \log(1 + \xi \frac{x_i}{\beta})}
#'
#' For all excess data points \eqn{n}:
#'
#' \deqn{l(\xi,\beta)=\sum_{i=1}^{n} (-\log(\beta) - (\frac{1}{\xi} + 1) \log(1 + \xi \frac{x_i}{\beta}))}
#'
#' \deqn{l(\xi,\beta)=-n\log(\beta) - (\frac{1}{\xi} + 1)\sum_{i=1}^{n}  \log(1 + \xi \frac{x_i}{\beta})}
#'
#' We can thus minimize \eqn{-l(\xi,\beta)}. The parameters \eqn{\xi} and \eqn{\beta} that minimize the negative log likelihood are the same that maximize the log likelihood. Hence, by using the excesses, we are able to determine \eqn{\xi} and \eqn{\beta} that best fit the tail of the data.
#'
#' There is also the case to consider when \eqn{\xi = 0} which results in an exponential distribution. The total log likelihood in such a case is:
#'
#' \deqn{l(0, \beta) =  -n \log(\beta) - \frac{1}{\beta} \sum_{i=1}^{n} x_i}
#'
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param u A numeric scalar that specifies the threshold value to calculate excesses
#' @param start_xi Initial value of \eqn{\xi} to pass to the optimizer
#' @param start_beta Initial value of \eqn{\beta} to pass to the optimizer
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns An unnamed numeric vector of length 2 containing the estimated
#' Generalized Pareto Distribution (GPD) parameters that minimize the negative log likelihood: \eqn{\xi} (shape/tail index) and \eqn{\beta} (scale parameter).
#' @export
#'
#' @examples
#' x <- rweibull(n=800, shape = 0.8, scale = 1)
#' values <- pot_estimator(data = x, u = 2, start_xi = 0.1, start_beta = NULL)
#'
#' @references
#'
#' Davison, A. C., & Smith, R. L. (1990). Models for exceedances over high thresholds. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{52}(3), 393-425. \doi{10.1111/j.2517-6161.1990.tb01796.x}
#'
#' Balkema, A. A., & de Haan, L. (1974). Residual life time at great age. \emph{The Annals of Probability}, \bold{2}(5), 792-804. \doi{10.1214/aop/1176996548}
#'
#' Pickands, J. (1975). Statistical Inference Using Extreme Order Statistics. \emph{The Annals of Statistics}, \bold{3}(1), 119–131. \url{http://www.jstor.org/stable/2958083}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 221-226) \doi{10.1017/9781009053730}
#'
pot_estimator <- function(data, u, start_xi=0.1, start_beta = NULL, na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(u) || length(u)!=1) {
    stop("`u` must be a numeric scalar")
  }

  if(u <= 0) {
    stop("Threshold u must be positive and within the support of the data")
  }

  if (!is.numeric(start_xi) || length(start_xi) != 1) {
    stop("`start_xi` must be a numeric scalar.")
  }
  if (!is.null(start_beta) && (!is.numeric(start_beta) || length(start_beta) != 1)) {
    stop("`start_beta` must be a numeric scalar.")
  }

  if (length(data) <= 1 ) {
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

  if (u >= max(data)) {
    stop("Threshold too high. No excesses possible.")
  }

  excesses <- data[data>u]- u

  if (length(excesses) < 2) {
    warning("Not enough data points above the threshold 'u' to fit the model.")
    return(c(xi = NA, beta = NA))
  }

  # Random starting points
  start_xi <- start_xi
  if (is.null(start_beta)) {
    start_beta <- sd(excesses, na.rm=TRUE)
  } else {
    start_beta <- start_beta
  }
  start_params <- c(xi = start_xi, beta = start_beta)

  opt_result <- optim(
    par = start_params,
    fn = gpd_lg_likelihood,
    data = excesses,
    method = "BFGS"
  )

  return(unname(opt_result$par)) #Removing names because plot() throws warnings otherwise
}
