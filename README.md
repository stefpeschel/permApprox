
# permApprox - Permutation p-value approximation

permApprox is an R package to compute empirical and approximated
p-values from permutation tests, especially useful when the number of
permutations is small and zero p-values should be strictly avoided. It
offers approximation using the Gamma distribution or the Generalized
Pareto Distribution (GPD) fitted to the tail of the permutation
distribution.

## Early version disclaimer

**Note:** This is an early version of the permApprox package and is
still under active development. While the core functionality is in
place, the interface, parameter settings, and behavior of certain
functions may change in future versions. Users are encouraged to use the
package with care and avoid relying on its current structure for
production workflows.

## Installation

``` r
install.packages("devtools")
devtools::install_github("stefpeschel/permApprox")
```

## Features

- Computes empirical p-values for permutation tests
- Approximates small p-values using the Gamma or GPD distributions
- Strictly avoids zero p-values
- Includes multiple testing correction support
- Control objects allow fine-tuning of each fitting method

## Usage

``` r
library(permApprox)

# Simulated example
set.seed(12345)
obs <- c(2.0, 3.0, 4.0, 5.0)
perm <- matrix(rnorm(4000), nrow = 4)

# Empirical p-values
res_emp <- compute_p_values(obs_stats = obs,
                            perm_stats = perm,
                            method = "empirical")

# Gamma approximation
gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
res_gamma <- compute_p_values(obs_stats = obs,
                              perm_stats = perm,
                              method = "gamma",
                              gamma_ctrl = gamma_ctrl)

# GPD approximation without constraints
gpd_ctrl <- make_gpd_ctrl(constraint = "unconstrained")
res_gpd <- compute_p_values(obs_stats = obs,
                            perm_stats = perm,
                            method = "gpd")

# GPD approximation with constraint
# (GPD must have support at t_obs + epsilon)
gpd_ctrl <- make_gpd_ctrl(constraint = "support_at_obs")
res_gpd_constr <- compute_p_values(obs_stats = obs,
                                   perm_stats = perm,
                                   method = "gpd",
                                   gpd_ctrl = gpd_ctrl)

# GPD approximation with constraint and fixed epsilon
gpd_ctrl <- make_gpd_ctrl(constraint = "support_at_max",
                          eps_fun = eps_fixed,
                          eps_par = list(value = 0.1))

res_gpd_constr_eps0.1 <- compute_p_values(obs_stats = obs,
                                   perm_stats = perm,
                                   method = "gpd",
                                   gpd_ctrl = gpd_ctrl)

# Data frame with (unadjusted) p-values
p_values <- data.frame(empirical = res_emp$p_unadjusted,
                       gamma = res_gamma$p_unadjusted,
                       gpd = res_gpd$p_unadjusted,
                       gpd_constr = res_gpd_constr$p_unadjusted,
                       gpd_constr_eps0.1 = res_gpd_constr_eps0.1$p_unadjusted)

p_values
```

    ##     empirical        gamma          gpd   gpd_constr gpd_constr_eps0.1
    ## 1 0.047952048 0.0596524811 4.316868e-02 4.441052e-02      4.317300e-02
    ## 2 0.005994006 0.0112178961 3.007131e-03 3.007131e-03      3.007131e-03
    ## 3 0.000999001 0.0023100226 2.416279e-05 1.489682e-08      2.319110e-05
    ## 4 0.000999001 0.0005513698 2.646884e-12 2.646884e-12      8.585641e-13

## Documentation

- The main function is `compute_p_values()`
- Control objects can be created using:
  - `make_gpd_ctrl()` for GPD fitting
  - `make_gamma_ctrl()` for Gamma fitting
  - `make_adjust_ctrl()` for multiple testing correction

Use `?compute_p_values` in R to see the full documentation and available
options.

## License

This package is licensed under GPL-3.
