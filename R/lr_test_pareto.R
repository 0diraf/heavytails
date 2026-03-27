#### Likelihood Ratio Test: Pareto vs. Alternatives ####

#' Likelihood Ratio Test: Pareto vs. Alternative Distributions
#'
#' Compares the Pareto distribution fit against one or more alternative
#' distributions using the Vuong likelihood ratio test for non-nested models
#' (§8.5, Step 3; Clauset et al. 2009, §3.3).
#'
#' For each alternative, the log-likelihood ratio
#' \eqn{LR = \ell_{\text{Pareto}} - \ell_{\text{alternative}}} is computed.
#' The Vuong test statistic checks whether the mean per-observation
#' log-likelihood ratio is significantly different from zero. A positive
#' \eqn{LR} with a small p-value indicates the Pareto is preferred; a negative
#' \eqn{LR} with a small p-value indicates the alternative is preferred.
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param xm A positive numeric scalar: lower bound. Only
#'   \code{data[data >= xm]} is used.
#' @param alternatives A character vector naming the distributions to compare
#'   against. Supported: \code{"exponential"}, \code{"lognormal"},
#'   \code{"weibull"}.
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#' @returns A \code{data.frame} with one row per alternative and columns:
#' \itemize{
#'   \item{\code{alternative}:} Name of the alternative distribution.
#'   \item{\code{ll_pareto}:} Pareto log-likelihood.
#'   \item{\code{ll_alternative}:} Alternative log-likelihood.
#'   \item{\code{lr_statistic}:} Vuong test statistic (z-score).
#'   \item{\code{p_value}:} Two-sided p-value.
#'   \item{\code{preferred}:} \code{"pareto"} or the alternative name.
#' }
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' x <- rpareto(n = 500, alpha = 2, xm = 1)
#' lr_test_pareto(x, xm = 1)
#'
#' @references
#'
#' Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions
#' in empirical data. \emph{SIAM Review}, \bold{51}(4), 661-703.
#' \doi{10.1137/070710111}
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, \bold{57}(2), 307-333.
#'
lr_test_pareto <- function(data, xm = NULL,
                           alternatives = c("exponential", "lognormal", "weibull"),
                           na.rm = FALSE) {

  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.null(xm) && (!is.numeric(xm) || length(xm) != 1 || xm <= 0)) {
    stop("`xm` must be a positive numeric scalar or NULL.")
  }

  valid_alts <- c("exponential", "lognormal", "weibull")
  if (!is.character(alternatives) || !all(alternatives %in% valid_alts)) {
    stop("`alternatives` must be a subset of: ", paste(valid_alts, collapse = ", "))
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

  x <- data[data >= xm]
  n <- length(x)

  if (n < 3) {
    stop("Not enough observations >= `xm` to run the test.")
  }

  # Fit Pareto via MLE
  pareto_fit <- mle_pareto(x, xm = xm, bias_corrected = FALSE)
  alpha_hat  <- pareto_fit$alpha

  # Per-observation Pareto log-likelihoods
  ll_obs_pareto <- dpareto(x, alpha = alpha_hat, xm = xm)
  ll_obs_pareto <- log(ll_obs_pareto[ll_obs_pareto > 0])
  ll_pareto <- sum(ll_obs_pareto)

  results <- lapply(alternatives, function(alt) {
    ll_obs_alt <- tryCatch(
      .fit_alternative_ll(x, alt),
      error = function(e) {
        warning(paste0("Failed to fit '", alt, "': ", e$message))
        rep(NA_real_, n)
      }
    )

    lr_obs <- ll_obs_pareto - ll_obs_alt
    lr_obs_valid <- lr_obs[is.finite(lr_obs)]
    n_valid <- length(lr_obs_valid)

    if (n_valid < 2) {
      return(data.frame(
        alternative = alt,
        ll_pareto = ll_pareto,
        ll_alternative = NA_real_,
        lr_statistic = NA_real_,
        p_value = NA_real_,
        preferred = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    # Vuong test statistic
    mean_lr <- mean(lr_obs_valid)
    sd_lr   <- stats::sd(lr_obs_valid)
    vuong_z <- sqrt(n_valid) * mean_lr / sd_lr
    p_val   <- 2 * stats::pnorm(-abs(vuong_z))

    preferred <- if (!is.finite(vuong_z)) NA_character_ else
      if (vuong_z > 0) "pareto" else alt

    data.frame(
      alternative = alt,
      ll_pareto = ll_pareto,
      ll_alternative = sum(ll_obs_alt[is.finite(ll_obs_alt)]),
      lr_statistic = vuong_z,
      p_value = p_val,
      preferred = preferred,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

# Internal: fit an alternative distribution and return per-obs log-likelihoods
#' @keywords internal
.fit_alternative_ll <- function(x, alt) {
  if (alt == "exponential") {
    rate_hat <- 1 / mean(x)
    stats::dexp(x, rate = rate_hat, log = TRUE)

  } else if (alt == "lognormal") {
    meanlog <- mean(log(x))
    sdlog   <- stats::sd(log(x))
    stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE)

  } else if (alt == "weibull") {
    # MLE for Weibull via optim on negative log-likelihood
    neg_ll_weibull <- function(params) {
      shape <- params[1]
      scale <- params[2]
      if (shape <= 0 || scale <= 0) return(1e10)
      -sum(stats::dweibull(x, shape = shape, scale = scale, log = TRUE))
    }
    # Starting values: method of moments approximation
    cv <- stats::sd(x) / mean(x)
    shape_init <- max(0.1, (cv)^(-1.086))
    scale_init <- mean(x) / gamma(1 + 1 / shape_init)
    opt <- stats::optim(par = c(shape_init, scale_init), fn = neg_ll_weibull,
                        method = "L-BFGS-B", lower = c(1e-6, 1e-6))
    stats::dweibull(x, shape = opt$par[1], scale = opt$par[2], log = TRUE)
  }
}
