#### Double Bootstrap ####

#' Estimators for Double Bootstrap.
#'
#' Helper function for doublebootstrap(). Calculates values for \eqn{\hat \xi^{(1)}}, and \eqn{M_{(n,k)}}.
#'
#' \deqn{\hat{\xi}_{n,k}^{(1)} := \frac{1}{k}\sum_{i=1}^{k}\log(\frac{X_{(i)}}{X_{(k+1)}})}
#'
#' \deqn{M_{n,k} := \frac{1}{k}\sum_{i=1}^{k} (\log(\frac{X_{(i)}}{X_{(k+1)}}))^2}
#'
#'
#' @param data Original data
#' @param k Number of top-order statistics
#'
#' @returns Returns a list containing values for \eqn{\hat \xi^{(1)}} determined through the Hill estimator, and \eqn{M_{(n,k)}}
#' @keywords internal
#'
#' @examples
#'
#' x <- heavytails:::generate_expbody_partail(n = 5000, alpha = 2.0)
#' db_est_values <- heavytails:::db_estimators(data = x, k = 10)
#'
db_estimators <- function(data, k) {

  x <- sort(data, decreasing = TRUE, na.last = NA)

  if (k <= 1 || k >= length(x) - 1) {
    return(list(hxi = NA, sec_moment = NA))
  }

  x_topk <- x[1:k]
  x_k1 <- x[k+1] #From Nair et al. (Eq. 9.12)

  if (x_k1 <= 0) {
    return(list(hxi = NA, sec_moment = NA))
  }

  log_ratios <- log(x_topk/x_k1)

  hxi <- mean(log_ratios)

  sec_moment <- mean(log_ratios[is.finite(log_ratios)]^2)

  return(list(hxi = hxi, sec_moment = sec_moment))
}

#' Difference calculation for Double Bootstrap
#'
#' Calculates the differences between \eqn{\hat \xi^{(1)}} and \eqn{\hat \xi^{(2)}} for each k value obtained through a sequence.
#'
#' Calculates \eqn{(M_{n,k}(j) - 2(\hat{\xi}_{n,k}^{(1)}(j))^2)^2}
#'
#' @param data Original dataset
#' @param k_seq Sequence of top-order statistics
#'
#' @returns Returns a vector of squared differences between \eqn{\hat \xi^{(1)}} and \eqn{M_{n,k}}
#' @keywords internal
#'
#' @examples
#'
#' x <- heavytails:::generate_expbody_partail(n = 5000, alpha = 2.0)
#' diffs <- heavytails:::find_db_diffs(data = x, k_seq = c(2,100))
#'
#'
find_db_diffs <- function(data, k_seq) {

  sq_diffs <- sapply(k_seq, function(k) {
    ests <- db_estimators(data, k)
    if (is.na(ests$hxi)) {
      return(NA)
    }

    diff <- ests$sec_moment-2*(ests$hxi^2)
    return(diff^2)
  })

  return(sq_diffs)
}
#' Double Bootstrap algorithm
#'
#' This function implements the Double Bootstrap algorithm as described by in Chapter 9 by *Nair et al.* It applies bootstrapping to two samples of different sizes to choose the value of \eqn{k} that minimizes the mean square error.
#'
#' Chapter 9 of *Nair et al.* specifically describes the Double Bootstrap algorithm for the Hill estimator.
#'
#' The Hill Double Bootstrap method uses the Hill estimator as the first estimator
#'
#' \deqn{\hat{\xi}_{n,k}^{(1)} := \frac{1}{k}\sum_{i=1}^{k}\log\left(\frac{X_{(i)}}{X_{(k+1)}}\right)}
#'
#' And a second moments-based estimator:
#'
#' \deqn{\hat{\xi}_{n,k}^{(2)} = \frac{M_{n,k}}{2\hat{\xi}_{n,k}^{H} }}
#'
#' Where
#'
#' \deqn{M_{n,k} := \frac{1}{k}\sum_{i=1}^{k}\left(\log\left(\frac{X_{(i)}}{X_{(k+1)}}\right)\right)^2}
#'
#' The difference between these two \eqn{\hat \xi} is given by:
#'
#' \deqn{|\hat{\xi}_{n,k}^{(1)} - \hat{\xi}_{n,k}^{(2)}| = \frac{|M_{n,k}-2(\hat{\xi}_{n,k}^{H})^{2}|}{2|\hat{\xi}_{n,k}^{H}|}}
#'
#' The Hill bootstrap method selects \eqn{\hat \kappa} in a way that minimizes the mean square error in the numerator by going through \eqn{r} bootstrap samples of different sizes \eqn{n_1} and \eqn{n_2}.
#'
#' \deqn{\hat{\kappa}_{1}^{*} := \text{arg min}_{k} \frac{1}{r} \sum_{j=1}^{r} (M_{n_1,k}(j) - 2(\hat{\xi}_{n_1,k}^{(1)}(j))^2)^2}
#'
#' This process is repeated to determine \eqn{\hat \kappa_{2}} with the bootstrap sample of size \eqn{n_{2}}. The final \eqn{\hat \kappa} is given by:
#'
#' \deqn{\hat{\kappa}^{*} = \frac{(\hat{\kappa}_{1}^{*})^2}{\hat{\kappa}_{2}^{*}} (\frac{\log \hat{\kappa}_{1}^{*}}{2\log n_1 - \log \hat{\kappa}_{1}^{*}})^{\frac{2(\log n_1 - \log \hat{\kappa}_{1}^{*})}{\log n_1}}}
#'
#' @param data A numeric vector of i.i.d. observations.
#' @param n1 A numeric scalar specifying the first bootstrap sample size, *Nair et al.* describe this as \eqn{n_1 = O(n^{1-\epsilon})} for \eqn{\epsilon \in (0, 1/2)}. Hence, default value (if n1 = -1) is chosen as 0.9.
#' @param n2 A numeric scalar specifying the second bootstrap sample size
#' @param r A numeric scalar specifying the number of bootstraps
#' @param k_max_prop A numeric scalar. The max k as a proportion of the sample size. It might be computationally expensive to consider all possible k values from the data. Furthermore, lower k values can be noisy, while higher values can be biased. Hence, k here is limited to a specific proportion (by default 50%) of the data
#' @param kvalues An integer specifying the length of sequence of candidate k values
#' @param na.rm Logical. If \code{TRUE}, missing values (\code{NA}) are removed
#'   before analysis. Defaults to \code{FALSE}.
#'
#'
#' @returns A named list containing the final results of the Double Bootstrap algorithm:
#' \itemize{
#'   \item{\code{k}:} The optimal number of top-order statistics \eqn{\hat{k}} selected by minimizing the MSE.
#'   \item{\code{alpha}:} The estimated tail index \eqn{\hat{\alpha}} (Hill estimator) corresponding to \eqn{\hat{k}}.
#' }
#'
#' @export
#'
#' @examples
#' xmin <- 1
#' alpha <- 2
#' r <- runif(800, 0, 1)
#' x <- (xmin * r^(-1/(alpha)))
#' db_kalpha <- doublebootstrap(data = x, n1 = -1, n2 = -1, r = 5, k_max_prop = 0.5, kvalues = 20)
#'
#' @references
#'
#' Danielsson, J., de Haan, L., Peng, L., & de Vries, C. G. (2001). Using a bootstrap method to choose the sample fraction in tail index estimation. \emph{Journal of Multivariate Analysis}, \bold{76}(2), 226–248. \doi{10.1006/jmva.2000.1903}
#'
#' Nair, J., Wierman, A., & Zwart, B. (2022). \emph{The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation}. Cambridge University Press. (pp. 229-233) \doi{10.1017/9781009053730}
#'
doublebootstrap <- function(data, n1 = -1, n2 = -1, r = 50, k_max_prop = 0.5, kvalues = 20, na.rm = FALSE) {


  if (!is.numeric(data) || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.")
  }

  if (!is.numeric(n1) || length(n1)!=1) {
    stop("`n1` must be a numeric scalar")
  }

  if (!is.numeric(n2) || length(n2)!=1) {
    stop("`n2` must be a numeric scalar.")
  }

  if (!is.numeric(r) || length(r)!=1) {
    stop("`r` must be numeric scalar.")
  }

  if (!is.numeric(k_max_prop) || length(k_max_prop) != 1 || k_max_prop <= 0 || k_max_prop >= 1) {
    stop("`k_max_prop` must be a numeric scalar in (0, 1).")
  }

  if (!is.numeric(kvalues) || kvalues < 5) {
    stop("`kvalues` must be numeric and >= 5.")
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

  if (all(data <= 0)) {
    stop("`data` must contain positive values for power-law fitting.")
  }

  if (any(data <= 0)) {
    warning("`data` has negative values. They will be excluded from calculations.")
    data <- data[data > 0]
    if (length(data) <= 5) {
      stop("`data` has length <= 5 after removing non-positive values.")
    }
  }


  n <- length(data)


  if (n1 == -1) {
    n1 <- floor(n^0.9)
  }

  if (n1 >= n)
    stop("'n1' must be less than the sample size n.")

  if (n1 <= 5) {
    stop("'n1' must be > 5.")
  }

  if (n2 == -1) {
    n2 <- floor(n1^2/n)
  }

  if (n2 >= n1)
    stop("'n2' must be less than n1.")

  if (n2 <= 5) {
    stop("'n2' must be > 5.")
  }

  k_max = floor(n1 * k_max_prop)

  if (k_max >= n1) {
    k_max <- n1 - 2
  }

  k_seq_n1 <- floor(seq(2, k_max, length.out = kvalues))


  all_sq_diffs_n1 <- matrix(NA, nrow = r, ncol = length(k_seq_n1))

  for (j in 1:r) {
    boot_sample_n1 <- sample(data, n1, replace = TRUE)
    all_sq_diffs_n1[j, ] <- find_db_diffs(boot_sample_n1, k_seq_n1)
  }

  mean_sq_diffs_n1 <- colMeans(all_sq_diffs_n1, na.rm = TRUE)

  k_hat_1 <- k_seq_n1[which.min(mean_sq_diffs_n1)]

  k_max_n2 <- floor(n2 * k_max_prop)
  if (k_max_n2 >= n2) {
    k_max_n2 <- n2 - 2
  }

  k_seq_n2 <- floor(seq(2, k_max_n2, length.out = kvalues))

  all_sq_diffs_n2 <- matrix(NA, nrow = r, ncol = length(k_seq_n2))

  for (j in 1:r) {
    boot_sample_n2 <- sample(data, n2, replace = TRUE)
    all_sq_diffs_n2[j, ] <- find_db_diffs(boot_sample_n2, k_seq_n2)
  }

  mean_sq_diffs_n2 <- colMeans(all_sq_diffs_n2, na.rm = TRUE)

  k_hat_2 <- k_seq_n2[which.min(mean_sq_diffs_n2)]

  k_hat_final <- floor(((k_hat_1^2) / k_hat_2) * (log(k_hat_1)/(2*log(n1)-log(k_hat_1)))^(2*(log(n1)-log(k_hat_1))/log(n1)))

  if(k_hat_final >= n) {
    k_hat_final <- n - 2
  }

  k_hat_final <- max(2, min(k_hat_final, n - 2))

  alpha_hat_final <- hill_estimator(data, k_hat_final) # To get the corresponding alpha_hat

  return(list(k=k_hat_final, alpha=alpha_hat_final))
}
