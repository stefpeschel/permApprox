
# permApprox <img src="man/figures/logo.png" align="right" height="180" />

The `permApprox` (**perm**utation p-value **approx**imation) package
provides tools to **approximate small permutation p-values** when the
number of permutations is limited and empirical p-values become coarse
or problematic (e.g., many ties or values very close to zero).

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

## Development version

Everyone who wants to use new features not included in any releases is
invited to install the development version:

``` r
remotes::install_github("stefpeschel/permApprox@develop")
```

Please check the
[NEWS](https://github.com/stefpeschel/permApprox/blob/develop/NEWS.md)
document for features implemented on develop branch.

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
X[group == 1, 1] <- X[group == 1, 1] + 1.5

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
    ## Number of tests             : 10
    ## Approximation method        : empirical (no approximation)
    ## Approximation threshold     : p-values < 0.1
    ## Multiple testing adjustment : BH
    ## 
    ## Final p-values:
    ##   min = 0.00999, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

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
    ## Number of tests             : 10
    ## Approximation method        : GPD tail approximation
    ## Approximation threshold     : p-values < 0.1
    ## Multiple testing adjustment : BH
    ## 
    ## Successful fits          : 3
    ## GOF rejections           : 0
    ## Fit failed               : 0
    ## No threshold found       : 0
    ## Discrete distributions   : 0
    ## Not selected for fitting : 7
    ## 
    ## Final p-values:
    ##   min = 0.0000000507, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

``` r
# More detailed diagnostics
summary(res_gpd)
```

    ## Summary of permApprox result
    ## ----------------------------
    ## Number of tests             : 10
    ## Approximation method        : GPD tail approximation
    ## Approximation threshold     : p-values < 0.1
    ## Multiple testing adjustment : BH
    ## 
    ## Fit status counts:
    ##   Successful fits          : 3
    ##   GOF rejections           : 0
    ##   Fit failed               : 0
    ##   No threshold found       : 0
    ##   Discrete distributions   : 0
    ##   Not selected for fitting : 7
    ## 
    ## GPD parameter summary (successful fits)
    ## --------------------------------------
    ##   shape:
    ##     min = -0.158, median = -0.12, mean = -0.12, max = -0.082
    ##   scale:
    ##     min = 0.113, median = 0.128, mean = 0.127, max = 0.141
    ##   n_exceed:
    ##     min =  250, median =  250, mean =  250, max =  250
    ## 
    ## P-value summary
    ## ---------------
    ## Empirical p-values:
    ##   empirical:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (unadjusted):
    ##   unadjusted:
    ##     min = 0.00000000507, median = 0.31, mean = 0.389, max = 0.987
    ## 
    ## Final p-values (adjusted, BH):
    ##   adjusted:
    ##     min = 0.0000000507, median = 0.562, mean = 0.548, max = 0.987
    ##   Rejections at alpha = 0.05: 1

``` r
# Inspect GPD fit results
fit_gpd <- res_gpd$fit_result

# Convert to a data frame
fit_gpd_df <- data.frame(
  pvals_emp     = res_gpd$p_empirical,
  pvals_gpd     = res_gpd$p_unadjusted,
  fit_status    = fit_gpd$status,
  shape         = fit_gpd$shape,
  thresh        = fit_gpd$thresh,
  n_exceed      = fit_gpd$n_exceed
)

fit_gpd_df
```

    ##      pvals_emp    pvals_gpd   fit_status       shape    thresh n_exceed
    ## 1  0.000999001 5.071830e-09      success -0.08196804 0.3180909      250
    ## 2  0.490509491 4.905095e-01 not_selected          NA        NA       NA
    ## 3  0.848151848 8.481518e-01 not_selected          NA        NA       NA
    ## 4  0.271728272 2.717283e-01 not_selected          NA        NA       NA
    ## 5  0.987012987 9.870130e-01 not_selected          NA        NA       NA
    ## 6  0.172827173 1.728272e-01 not_selected          NA        NA       NA
    ## 7  0.348651349 3.486513e-01 not_selected          NA        NA       NA
    ## 8  0.652347652 6.523477e-01 not_selected          NA        NA       NA
    ## 9  0.040959041 3.981341e-02      success -0.15785650 0.2496041      250
    ## 10 0.081918082 8.289617e-02      success -0.11966759 0.2653776      250

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
    ## Number of tests             : 10
    ## Approximation method        : Gamma approximation
    ## Approximation threshold     : p-values < 0.1
    ## Multiple testing adjustment : BH
    ## 
    ## Successful fits          : 3
    ## GOF rejections           : 0
    ## Fit failed               : 0
    ## Discrete distributions   : 0
    ## Not selected for fitting : 7
    ## 
    ## Final p-values:
    ##   min = 0.000862, median = 0.562, max = 0.987
    ## 
    ## Use summary() for detailed fit diagnostics.

``` r
summary(res_gamma)
```

    ## Summary of permApprox result
    ## ----------------------------
    ## Number of tests             : 10
    ## Approximation method        : Gamma approximation
    ## Approximation threshold     : p-values < 0.1
    ## Multiple testing adjustment : BH
    ## 
    ## Fit status counts:
    ##   Successful fits          : 3
    ##   GOF rejections           : 0
    ##   Fit failed               : 0
    ##   Discrete distributions   : 0
    ##   Not selected for fitting : 7
    ## 
    ## Gamma parameter summary (successful fits)
    ## ----------------------------------------
    ##   shape:
    ##     min = 1.33, median = 1.39, mean = 1.37, max =  1.4
    ##   rate:
    ##     min = 6.39, median = 7.21, mean = 7.25, max = 8.16
    ## 
    ## P-value summary
    ## ---------------
    ## Empirical p-values:
    ##   empirical:
    ##     min = 0.000999, median = 0.31, mean = 0.39, max = 0.987
    ## 
    ## Final p-values (unadjusted):
    ##   unadjusted:
    ##     min = 0.0000862, median = 0.31, mean = 0.393, max = 0.987
    ## 
    ## Final p-values (adjusted, BH):
    ##   adjusted:
    ##     min = 0.000862, median = 0.562, mean = 0.563, max = 0.987
    ##   Rejections at alpha = 0.05: 1

``` r
# Inspect Gamma fit results
fit_gamma <- res_gamma$fit_result

# Convert to a data frame
fit_gamma_df <- data.frame(
  pvals_emp       = res_gamma$p_empirical,
  pvals_gamma     = res_gamma$p_unadjusted,
  fit_status      = fit_gamma$status,
  shape           = fit_gamma$shape,
  rate            = fit_gamma$rate
)

fit_gamma_df
```

    ##      pvals_emp  pvals_gamma   fit_status    shape     rate
    ## 1  0.000999001 8.616786e-05      success 1.403064 6.392932
    ## 2  0.490509491 4.905095e-01 not_selected       NA       NA
    ## 3  0.848151848 8.481518e-01 not_selected       NA       NA
    ## 4  0.271728272 2.717283e-01 not_selected       NA       NA
    ## 5  0.987012987 9.870130e-01 not_selected       NA       NA
    ## 6  0.172827173 1.728272e-01 not_selected       NA       NA
    ## 7  0.348651349 3.486513e-01 not_selected       NA       NA
    ## 8  0.652347652 6.523477e-01 not_selected       NA       NA
    ## 9  0.040959041 6.044530e-02      success 1.390864 8.162339
    ## 10 0.081918082 9.806895e-02      success 1.325764 7.207294

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
    ## 1    5.071830e-09 3.590799e-08
    ## 2    4.905095e-01 4.961075e-01
    ## 3    8.481518e-01 6.672022e-01
    ## 4    2.717283e-01 3.847611e-01
    ## 5    9.870130e-01 6.987941e-01
    ## 6    1.728272e-01 3.058993e-01
    ## 7    3.486513e-01 4.114021e-01
    ## 8    6.523477e-01 5.773185e-01
    ## 9    3.981341e-02 1.409372e-01
    ## 10   8.289617e-02 1.956319e-01

------------------------------------------------------------------------

## License

This package is released under the GPL-3 license.
