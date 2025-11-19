
# permApprox <img src="man/figures/logo.png" align="right" height="180" />

`permApprox` provides tools to **approximate small permutation
p-values** when the number of permutations is limited and empirical
p-values become coarse or problematic (e.g., many ties or values very
close to zero).

The core function `perm_approx()` computes empirical p-values and, for
tests with small empirical p-values, replaces them by parametric
approximations derived from the permutation statistics using either

- a **Generalized Pareto Distribution (GPD)** fitted to the **upper
  tail** (exceedances above a data-driven threshold), or  
- a **Gamma distribution** fitted to the **full permutation
  distribution**.

For the GPD approach, `permApprox` supports **support constraints** to
avoid zero p-values when the estimated shape parameter is negative.
Multiple testing adjustment is performed on the resulting p-values using
standard p-value based methods (e.g. BH, adaptive BH, local FDR).

------------------------------------------------------------------------

## Installation

The package is currently available from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("stefpeschel/permApprox")
```

Load the package:

``` r
library(permApprox)
```

------------------------------------------------------------------------

## Overview

Permutation tests are widely used when the null distribution of a test
statistic is not available in closed form. With a limited number of
permutations, empirical p-values are often too coarse and may even be
zero for strong effects, which causes problems for multiple testing
procedures.

`permApprox` addresses this by:

1.  Computing **empirical p-values** from permutation test statistics,
2.  Selecting tests with empirical p-values below an **approximation
    threshold**,
3.  Fitting a **GPD** to the tail (or a **Gamma** distribution to the
    full permutation distribution),
4.  Enforcing **support constraints** to avoid zero p-values for
    negative shape parameters, and
5.  Optionally applying **multiple testing correction** (BH, adaptive
    BH, lfdr, …).

The main entry point is:

``` r
perm_approx(
  obs_stats,    # observed test statistics
  perm_stats,   # permutation test statistics (rows = permutations, cols = tests)
  method = c("gpd", "gamma", "empirical"),
  ...
)
```

------------------------------------------------------------------------

## Usage

### Basic example: 10 tests, one true effect

We illustrate the workflow with a simple two-sample mean comparison for
10 tests, where the first test has a true effect.

``` r
set.seed(42)

n_per_group <- 50
m_tests     <- 10

# Group labels: 0 = control, 1 = treatment
group <- rep(c(0, 1), each = n_per_group)

# Data matrix: rows = samples, cols = tests/features
X <- matrix(
  rnorm(2 * n_per_group * m_tests),
  ncol = m_tests
)

# Introduce a true effect in the first test
X[group == 1, 1] <- X[group == 1, 1] + 0.8

# Observed test statistics: mean difference (treated - control)
obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
             colMeans(X[group == 0, , drop = FALSE])

# Permutation distribution: shuffle group labels and recompute all 10 stats
B <- 1000
perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)

for (b in seq_len(B)) {
  grp_perm <- sample(group)
  perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
                   colMeans(X[grp_perm == 0, , drop = FALSE])
}
```

### Empirical p-values

``` r
res_emp <- perm_approx(
  obs_stats  = obs_stats,
  perm_stats = perm_mat,
  method     = "empirical",
  verbose    = FALSE
)

res_emp
```

    ## permApprox result
    ## -----------------
    ## Number of tests              : 10
    ## Approximation method         : empirical (no approximation)
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Final p-values:
    ##   min = 0.00999, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

``` r
summary(res_emp)
```

    ## Summary of permApprox result
    ## ----------------------------
    ## Number of tests              : 10
    ## Approximation method         : empirical (no parametric tail approximation)
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Fit status counts:
    ##   (no parametric fit information available)
    ## 
    ## P-value summary
    ## ---------------
    ## Empirical p-values:
    ##   empirical:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (unadjusted):
    ##   unadjusted:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (adjusted, BH):
    ##   adjusted:
    ##     min = 0.00999, median = 0.562, mean = 0.549, max = 0.987
    ##   Rejections at alpha = 0.05: 1

``` r
# Adjusted (BH) and unadjusted p-values
res_emp$p_values
```

    ##  [1] 0.00999001 0.70072784 0.94239094 0.54345654 0.98701299 0.43206793
    ##  [7] 0.58108558 0.81543457 0.20479520 0.27306027

``` r
res_emp$p_unadjusted
```

    ##  [1] 0.000999001 0.490509491 0.848151848 0.271728272 0.987012987 0.172827173
    ##  [7] 0.348651349 0.652347652 0.040959041 0.081918082

### GPD approximation

We now approximate small p-values using a GPD tail model with a support
constraint at the observed statistic (single-test mode).

``` r
gpd_ctrl <- make_gpd_ctrl(
  constraint   = "support_at_obs",
  sample_size  = n_per_group
)

res_gpd <- perm_approx(
  obs_stats  = obs_stats,
  perm_stats = perm_mat,
  method     = "gpd",
  gpd_ctrl   = gpd_ctrl,
  verbose    = FALSE
)

# Compact overview
res_gpd
```

    ## permApprox result
    ## -----------------
    ## Number of tests              : 10
    ## Approximation method         : GPD tail approximation
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Successful fits              : 3
    ## Discrete distributions       : 0
    ## GOF rejections               : 0
    ## 
    ## Final p-values:
    ##   min = 0.000726, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

``` r
# More detailed diagnostics
summary(res_gpd)
```

    ## Summary of permApprox result
    ## ----------------------------
    ## Number of tests              : 10
    ## Approximation method         : GPD tail approximation
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Fit status counts:
    ##   success      : 3
    ##   gof_reject   : 0
    ##   fit_failed   : 0
    ##   discrete     : 0
    ##   no_threshold : 0
    ##   not_selected : 7
    ## 
    ## GPD parameter summary (successful fits)
    ## --------------------------------------
    ##   shape:
    ##     min = -0.158, median = -0.122, mean = -0.133, max = -0.12
    ##   scale:
    ##     min = 0.113, median = 0.128, mean = 0.124, max = 0.129
    ##   n_exceed:
    ##     min =  250, median =  250, mean =  250, max =  250
    ## 
    ## Goodness-of-fit p-values (all fitted tests)
    ## ------------------------------------------
    ##   GOF p-values:
    ##     min = 0.218, median = 0.533, mean = 0.551, max = 0.902
    ## 
    ## P-value summary
    ## ---------------
    ## Empirical p-values:
    ##   empirical:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (unadjusted):
    ##   unadjusted:
    ##     min = 0.0000726, median = 0.31, mean = 0.389, max = 0.987
    ## 
    ## Final p-values (adjusted, BH):
    ##   adjusted:
    ##     min = 0.000726, median = 0.562, mean = 0.548, max = 0.987
    ##   Rejections at alpha = 0.05: 1

``` r
# Tail fit details (from .compute_pvals_gpd())
fit_gpd <- res_gpd$fit_result
fit_gpd$status      # success / discrete / no_threshold / gof_reject / fit_failed
```

    ##  [1] success      not_selected not_selected not_selected not_selected
    ##  [6] not_selected not_selected not_selected success      success     
    ## Levels: not_selected discrete no_threshold fit_failed gof_reject success

``` r
fit_gpd$shape       # GPD shape parameters
```

    ##  [1] -0.1220435         NA         NA         NA         NA         NA
    ##  [7]         NA         NA -0.1578565 -0.1196676

``` r
fit_gpd$thresh      # thresholds per test
```

    ##  [1] 0.2688005        NA        NA        NA        NA        NA        NA
    ##  [8]        NA 0.2496041 0.2653776

``` r
fit_gpd$n_exceed    # number of exceedances
```

    ##  [1] 250  NA  NA  NA  NA  NA  NA  NA 250 250

Compare empirical and GPD-based p-values:

``` r
data.frame(
  empirical = res_emp$p_unadjusted,
  GPD       = res_gpd$p_unadjusted
)
```

    ##      empirical          GPD
    ## 1  0.000999001 0.0000726215
    ## 2  0.490509491 0.4905094905
    ## 3  0.848151848 0.8481518482
    ## 4  0.271728272 0.2717282717
    ## 5  0.987012987 0.9870129870
    ## 6  0.172827173 0.1728271728
    ## 7  0.348651349 0.3486513487
    ## 8  0.652347652 0.6523476523
    ## 9  0.040959041 0.0398134115
    ## 10 0.081918082 0.0828961722

### Gamma approximation

Alternatively, you can fit a Gamma distribution to the permutation
statistics:

``` r
gamma_ctrl <- make_gamma_ctrl(gof_test = "none")

res_gamma <- perm_approx(
  obs_stats   = obs_stats,
  perm_stats  = perm_mat,
  method      = "gamma",
  gamma_ctrl  = gamma_ctrl,
  verbose     = FALSE
)

res_gamma
```

    ## permApprox result
    ## -----------------
    ## Number of tests              : 10
    ## Approximation method         : Gamma approximation
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Successful fits              : 3
    ## Discrete distributions       : 0
    ## GOF rejections               : 0
    ## 
    ## Final p-values:
    ##   min = 0.023, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

``` r
summary(res_gamma)
```

    ## Summary of permApprox result
    ## ----------------------------
    ## Number of tests              : 10
    ## Approximation method         : Gamma approximation
    ## Approximation threshold      : p-values < 0.1
    ## Multiple testing adjustment  : BH
    ## 
    ## Fit status counts:
    ##   success      : 3
    ##   gof_reject   : 0
    ##   fit_failed   : 0
    ##   discrete     : 0
    ##   not_selected : 7
    ## 
    ## Gamma parameter summary (successful fits)
    ## ----------------------------------------
    ##   shape:
    ##     min = 1.33, median = 1.39, mean = 1.38, max = 1.43
    ##   rate:
    ##     min = 7.21, median = 7.58, mean = 7.65, max = 8.16
    ## 
    ## P-value summary
    ## ---------------
    ## Empirical p-values:
    ##   empirical:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (unadjusted):
    ##   unadjusted:
    ##     min = 0.0023, median = 0.31, mean = 0.393, max = 0.987
    ## 
    ## Final p-values (adjusted, BH):
    ##   adjusted:
    ##     min = 0.023, median = 0.562, mean = 0.565, max = 0.987
    ##   Rejections at alpha = 0.05: 1

``` r
gamma_fit <- res_gamma$fit_result
gamma_fit$shape
```

    ##  [1] 1.431798       NA       NA       NA       NA       NA       NA       NA
    ##  [9] 1.390864 1.325764

``` r
gamma_fit$rate
```

    ##  [1] 7.580890       NA       NA       NA       NA       NA       NA       NA
    ##  [9] 8.162339 7.207294

### Multiple testing adjustment

`permApprox` can pass empirical / approximated p-values to a
resampling-based multiple testing procedure. Supported methods include:

- `"none"` – no adjustment
- classical `stats::p.adjust()` methods (e.g. `"BH"`, `"holm"`, `"BY"`,
  …)
- `"lfdr"` – local FDR via `fdrtool`
- `"adapt_BH"` – adaptive BH based on estimated proportion of true nulls

Example using adaptive BH on GPD-based p-values:

``` r
adjust_ctrl <- make_adjust_ctrl(true_null_method = "lfdr")

res_adapt <- perm_approx(
  obs_stats     = obs_stats,
  perm_stats    = perm_mat,
  method        = "gpd",
  gpd_ctrl      = gpd_ctrl,
  adjust_method = "adapt_BH",
  adjust_ctrl   = adjust_ctrl,
  verbose       = FALSE
)

data.frame(
  GPD_unadjusted = res_gpd$p_unadjusted,
  GPD_adjusted   = res_adapt$p_values
)
```

    ##    GPD_unadjusted GPD_adjusted
    ## 1    0.0000726215 0.0005141617
    ## 2    0.4905094905 0.4961167117
    ## 3    0.8481518482 0.6672146680
    ## 4    0.2717282717 0.3847683174
    ## 5    0.9870129870 0.6988071646
    ## 6    0.1728271728 0.3059049582
    ## 7    0.3486513487 0.4114097511
    ## 8    0.6523476523 0.5773292997
    ## 9    0.0398134115 0.1409398741
    ## 10   0.0828961722 0.1956355179

------------------------------------------------------------------------

## Key features

- **Unified interface** for empirical, GPD, and Gamma-based p-values via
  `perm_approx()`.

- **Support-constrained GPD fitting** to avoid zero p-values for
  negative shape parameters, with options:

  - `constraint = "support_at_obs"` – per-test constraint at the
    observed statistic
  - `constraint = "support_at_max"` – global constraint at the maximum
    test statistic
  - `constraint = "unconstrained"`

- **Flexible epsilon definition** (evaluation point above the boundary)
  via user-specified epsilon functions.

- **Discreteness screening** for permutation distributions that are too
  coarse for reliable tail fitting.

- **Goodness-of-fit testing** (e.g. Cramér–von Mises) with fallback to
  empirical p-values if the tail model is rejected.

------------------------------------------------------------------------

## License

This package is released under the GPL-3 license.
