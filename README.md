# heavytails

**Estimators, diagnostics, and goodness-of-fit tools for heavy-tailed distributions in R.**

`heavytails` implements the estimators and algorithms from *The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation* (Nair, Wierman & Zwart, 2022) — Chapters 8 and 9. It covers tail index estimation, visual diagnostics, Pareto model fitting, and a complete goodness-of-fit pipeline, all using base R with no heavy dependencies.

## Installation

From CRAN:

```r
install.packages("heavytails")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("0diraf/heavytails")
```

## Quick Example

```r
library(heavytails)

set.seed(1)
x <- rpareto(n = 1000, alpha = 2, xm = 1)

# Tail index estimation
hill_estimator(x, k = 50)
mle_pareto(x)

# Visual diagnostics
hill_plot(x)
rank_plot(x)

# Goodness-of-fit pipeline
fit  <- plfit(x)
xmin <- fit$xmin
alph <- fit$alpha

ks_gof(x, alpha = alph, xm = xmin, n_boot = 500)
lr_test_pareto(x, alpha = alph, xmin = xmin)
```

## Functions

### Tail Index Estimators

| Function | Description | Reference |
|---|---|---|
| `hill_estimator()` | Hill (1975) MLE on top-k order statistics | Ch. 9, §9.2 |
| `moments_estimator()` | Dekkers-Einmahl-de Haan (1989) moments estimator | Ch. 9, §9.3 |
| `pickands_estimator()` | Pickands (1975) spacing-based estimator | Ch. 9, §9.3 |
| `pot_estimator()` | GPD fit via MLE on excesses over threshold u | Ch. 9, §9.4 |
| `plfit()` | KS-minimization for optimal k̂ and α̂ (Clauset et al. 2009) | Ch. 9, §9.5 |
| `doublebootstrap()` | Automatic k selection via double bootstrap (Danielsson et al. 2001) | Ch. 9, §9.5 |

### Pareto Estimation

| Function | Description |
|---|---|
| `mle_pareto()` | Parametric MLE for Pareto(xm, α) with optional bias correction |
| `wls_pareto()` | Weighted least-squares log-rank regression estimator |
| `ks_xmin()` | KS-based selection of the optimal xmin threshold |

### Visual Diagnostics

| Function | Description |
|---|---|
| `hill_plot()` | Hill plot: α̂ vs. k to assess stability |
| `moments_plot()` | Moments estimator plot: ξ̂ vs. k |
| `pickands_plot()` | Pickands estimator plot: ξ̂ vs. k |
| `rank_plot()` | Log-rank vs. log-x plot for Pareto linearity |
| `qq_pareto()` | Pareto Q-Q plot |

All plot functions return the underlying `data.frame` invisibly, so results can be captured for custom plotting with ggplot2 or base R.

### Goodness-of-Fit

| Function | Description |
|---|---|
| `ks_gof()` | Bootstrap KS test for Pareto fit |
| `lr_test_pareto()` | Likelihood-ratio test: Pareto vs. alternative distributions |

### Pareto Utilities

| Function | Description |
|---|---|
| `rpareto()` | Generate Pareto(xm, α) random variates |
| `pareto_cdf()` | Pareto CDF |
| `dpareto()` | Pareto density |

## Example: Hill Plot

A Hill plot shows how the tail index estimate changes with k. Stability across a range of k values is evidence that the tail is power-law distributed.

```r
library(heavytails)

set.seed(42)
x <- rpareto(n = 2000, alpha = 1.5, xm = 1)

hill_plot(x, alpha_true = 1.5,
          main = "Hill Plot", col = "steelblue")
```

The dashed red line (when `alpha_true` is supplied) marks the true value for simulation studies.

## Example: Full Clauset et al. Pipeline

```r
library(heavytails)

set.seed(1)
x <- rpareto(n = 1000, alpha = 2, xm = 1)

# Step 1: estimate xmin and alpha
fit <- plfit(x)

# Step 2: goodness-of-fit test
gof <- ks_gof(x, alpha = fit$alpha, xm = fit$xmin, n_boot = 1000)
gof$p_value   # > 0.1 → cannot reject Pareto

# Step 3: compare against alternative distributions
lr_test_pareto(x, alpha = fit$alpha, xmin = fit$xmin)
```

## Reference

Nair, J., Wierman, A., & Zwart, B. (2022). *The Fundamentals of Heavy Tails: Properties, Emergence, and Estimation*. Cambridge University Press. [doi:10.1017/9781009053730](https://doi.org/10.1017/9781009053730)

## License

MIT
