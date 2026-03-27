# heavytails 0.2.0

## New Functions

### Pareto Estimation

- `mle_pareto()` — Parametric MLE for the Pareto tail index using the closed-form estimator from Theorem 8.1 of Nair et al. Supports an optional finite-sample bias correction (§8.3). Returns a `heavytails_mle` object with a `print()` method.
- `wls_pareto()` — Weighted least-squares estimator for the Pareto tail index via log-log rank regression (Theorem 8.5 of Nair et al.). Optionally draws the rank plot with both WLS and OLS fitted lines. Also returns the unweighted OLS estimate for comparison.
- `ks_xmin()` — Standalone KS-minimization estimator for the power-law lower bound x̂_m (Step 1 of the Clauset et al. pipeline, §8.5). Extracts the core loop from `plfit()` so x̂_m can be estimated without running the full power-law fit.

### Goodness-of-Fit

- `ks_gof()` — Bootstrap Kolmogorov-Smirnov goodness-of-fit test for the Pareto distribution (Step 2 of the Clauset et al. pipeline, §8.5). The p-value is the fraction of parametric bootstrap KS statistics exceeding the observed statistic.
- `lr_test_pareto()` — Vuong likelihood-ratio test comparing the Pareto fit against alternative distributions: `"exponential"`, `"lognormal"`, and `"weibull"` (Step 3 of the Clauset et al. pipeline; Clauset et al. 2009, §3.3). Returns a `data.frame` with the LR statistic, two-sided p-value, and preferred model for each alternative.

### Visual Diagnostics

- `hill_plot()` — Hill plot: α̂ vs. k over a range of tail sizes. Accepts an optional `alpha_true` argument to overlay the true value as a reference line (useful in simulation studies). Returns the plotted `data.frame` invisibly.
- `moments_plot()` — Moments estimator plot: ξ̂ vs. k. Returns the plotted `data.frame` invisibly.
- `pickands_plot()` — Pickands estimator plot: ξ̂ vs. k. Returns the plotted `data.frame` invisibly.
- `rank_plot()` — Log-rank vs. log-x plot for assessing Pareto linearity. Returns the plotted `data.frame` invisibly.
- `qq_pareto()` — Pareto Q-Q plot. Accepts an optional `alpha` argument; if omitted, `mle_pareto()` is called internally on the data. Returns the theoretical and empirical quantiles invisibly.

### Pareto Utilities (newly exported)

- `rpareto()` — Generate Pareto(xm, α) random variates. Promoted from internal helper to exported function.
- `pareto_cdf()` — Pareto CDF. Promoted from internal helper to exported function.
- `dpareto()` — Pareto density function. New function.

## Bug Fixes

- **Hill estimator indexing corrected.** The denominator in the Hill formula was `X_(k)` instead of `X_(k+1)`, causing the k-th term in the average to always be zero and biasing the estimate upward. The formula now correctly divides by `X_(k+1)`, consistent with Eq. 9.12 of Nair et al. and the internal `db_estimators()` helper used by `doublebootstrap()`.
- **Dead code removed from `doublebootstrap()`.** The `if (n2 < 5)` check was unreachable because the preceding `if (n2 <= 5)` already covered that case. The `if (n1 != -1 && n1 >= n)` check was also unreachable because `n1` had already been overwritten from its `-1` default earlier in the function body. Both branches have been removed.

## Breaking Changes

- `gpd_lg_likelihood()` is no longer exported. It is an internal optimization helper for `pot_estimator()` and was never intended for direct use. Any code calling `gpd_lg_likelihood()` directly should switch to `heavytails:::gpd_lg_likelihood()` or, preferably, use `pot_estimator()` instead.

## Documentation

- Added `README.md` with installation instructions, a quick-start example, a function reference table, and a full Clauset et al. pipeline walkthrough.
- All new functions have full roxygen2 documentation with LaTeX math, parameter descriptions, return value specifications, examples, and references to Nair et al. (2022) and Clauset et al. (2009).
- Added `NEWS.md` (this file).

---

# heavytails 0.1.0

- Initial CRAN release.
- Implemented tail index estimators: `hill_estimator()`, `moments_estimator()`, `pickands_estimator()`, `pot_estimator()`, `plfit()`, `doublebootstrap()`.
