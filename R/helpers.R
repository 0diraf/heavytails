#### Helper Functions for Testing ####

#' @keywords internal
test_est <- function(data=NULL, seed=1, n=2000, alpha=2, xmin=1, kmax=2000, estimator, wshape=0.8, wscale=1) {

  set.seed(seed)
  if (is.null(data)) {      #This if check is for testing PLFIT and Double-bootstrap
    x <- rpareto(n = n, alpha = alpha, xm=xmin)
  } else {
    x <- data
  }


  if (estimator=='hill'){
    values <- sapply(1:kmax, function(k) hill_estimator(x, k))

  } else if (estimator=='moments') {
    values <- sapply(1:kmax, function(k) moments_estimator(x, k))
  } else if (estimator=='pickands') {
    values <- sapply(1:kmax, function(k) pickands_estimator(x, k))
  } else if (estimator=='pot') {
    y <- rweibull(n, shape = wshape, scale = wscale)
    u_max <- quantile(y, 0.999)
    u_values <- seq(0.1, u_max, length.out = 100)
    values <- sapply(u_values, function(u) pot_estimator(y, u))
    l <- list(th=u_values, v=values)
    return(l)
  }

  return(values)
}

#' @keywords internal
rexp_truncated <- function(n, rate, cutoff) {
  if (n == 0) {
    return(numeric(0))
  }
  samples <- numeric(n)
  count <- 0
  while (count < n) {

    batch_size <- 2*(n-count)
    candidates <- rexp(batch_size, rate = rate)

    good_samples <- candidates[candidates <= cutoff]

    num_to_add <- min(length(good_samples), n-count)

    if (num_to_add > 0) {
      samples[(count + 1):(count + num_to_add)] <- good_samples[1:num_to_add]
      count <- count + num_to_add
    }
  }

  return(samples)
}

#' @keywords internal
generate_expbody_partail <- function(n, prob_tail = 0.5,
                                     cutoff = 10, alpha = 2.0, rate = 0.2) {

  n_tail <- rbinom(1, size=n, prob=prob_tail)
  n_body <- n - n_tail

  data_tail <- rpareto(n =n_tail, xm=cutoff, alpha=alpha)
  data_body <- rexp_truncated(n=n_body, rate=rate, cutoff=cutoff)

  combined_data <- c(data_tail, data_body)

  return(sample(combined_data))
}
