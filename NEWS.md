
## permApprox 1.0.1 <img src="man/figures/logo.png" align="right" height="180" />

### 🐛 Bug-fixes

- `.find_gpd_thresh()` threw an error (“Threshold method not
  supported.”) if some of the GoF tests were accepted, but never three
  in a row. In these cases, the loop does not stop early.

## permApprox 1.0.0

### 🆕 New features

- Renamed the main function from `compute_p_values()` to
  **`perm_approx()`** — more concise and package-branded.
- Added a new argument **`approx_thresh`** (formerly `fit_thresh`) to
  specify the empirical p-value threshold under which parametric
  approximation is applied.
- Introduced **`print.perm_approx()`** and **`summary.perm_approx()`**
  S3 methods for cleaner interactive output and diagnostic summaries.
- Unified approximation output under a single slot **`fit_result`**;
  removed dual `gpd_fit` / `gamma_fit` slots for simplicity.
- Added a **`method_used`** output vector to track whether each test
  used empirical, Gamma, or GPD p-value.
- Extended the Gamma-approximation branch of the code: now includes a
  `status` factor (success, discrete, fit_failed, gof_reject) and
  parallel fitting support.
- Enhanced the GPD branch: improved threshold-detection methods, new
  `eps_factor()` handling, improved support-constraint handling
  (`support_at_obs`, `support_at_max`, `unconstrained`), and fixed
  edge-cases (shape = 0, zero exceedances).
- Added `cores` and `parallel_min` arguments in the main function for
  parallel execution of threshold detection and fitting (via the
  `future` framework).
- Cleaned up and simplified the workflow: one-sided transformation moved
  to global helper, argument defaults updated, improved error-handling
  and documentation.
- Updated examples in documentation to reflect the new API (including
  print/summary calls for users).

### 🐛 Bug-fixes

- Fixed internal `\_pgpd_upper_tail()` behavior for shape = 0.
- Fixed the `.compute_pvals_gamma()` logic to correctly calculate scaled
  tail probabilities when `include_obs = TRUE`.
- Changed default parameter values in control objects to more sensible
  settings; fixed mismatches between documentation and code.
- Several improvements to `threshold_detection` logic: cleaned up method
  names, improved stability of exceedance counts, and fixed progress-bar
  behaviour on Windows.

### 📚 Documentation & tests

- Updated the README, manual pages, and examples to reflect the new
  features.
- Added comprehensive `testthat` tests covering empirical mode, GPD
  approximation, Gamma approximation, print/summary methods, and
  multiple-testing adjustment behaviour.
- Streamlined package imports and removed unused dependencies.
