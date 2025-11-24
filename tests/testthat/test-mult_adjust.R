# tests/testthat/test-mult_adjust.R

set.seed(42)

## Helper: simple simulation of m tests with some signal
simulate_multi <- function(m = 50, n = 100, B = 500, n_signal = 3, effect = 1.5) {
  group <- rep(c(0, 1), each = n)
  
  X <- matrix(rnorm(2 * n * m), ncol = m)
  if (n_signal > 0) {
    for (j in seq_len(n_signal)) {
      X[group == 1, j] <- X[group == 1, j] + effect
    }
  }
  
  obs <- colMeans(X[group == 1, ]) - colMeans(X[group == 0, ])
  
  perm <- matrix(NA_real_, nrow = B, ncol = m)
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm[b, ] <- colMeans(X[grp_perm == 1, ]) - colMeans(X[grp_perm == 0, ])
  }
  
  list(obs = obs, perm = perm)
}

## ---------------------------------------------------------------------------
## 1. adjust_method = "none"
## ---------------------------------------------------------------------------

test_that("adjust_method = 'none' returns unadjusted p-values", {
  sim <- simulate_multi()
  res <- perm_approx(
    obs_stats     = sim$obs,
    perm_stats    = sim$perm,
    method        = "empirical",
    adjust_method = "none",
    verbose       = FALSE
  )
  
  expect_true(all(res$p_values == res$p_unadjusted))
  expect_null(res$adjust_result)
})

## ---------------------------------------------------------------------------
## 2. Classical p.adjust() methods: BH, holm, BY
## ---------------------------------------------------------------------------

test_that("Classical p.adjust methods work (BH, holm, BY)", {
  sim <- simulate_multi()
  methods <- c("BH", "holm", "BY")
  
  for (meth in methods) {
    res <- perm_approx(
      obs_stats     = sim$obs,
      perm_stats    = sim$perm,
      method        = "empirical",
      adjust_method = meth,
      verbose       = FALSE
    )
    
    expect_s3_class(res, "perm_approx")
    expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
    expect_true(all(res$p_values >= res$p_unadjusted))
    expect_false(is.null(res$adjust_result))
  }
})

## ---------------------------------------------------------------------------
## 3. lfdr-based adjustment
## ---------------------------------------------------------------------------

test_that("lfdr adjustment works and returns valid structure", {
  sim <- simulate_multi(m = 200, n = 100, n_signal = 5)
  
  res <- perm_approx(
    obs_stats     = sim$obs,
    perm_stats    = sim$perm,
    method        = "empirical",
    adjust_method = "lfdr",
    adjust_ctrl   = make_adjust_ctrl(true_null_method = "lfdr"),
    verbose       = FALSE
  )
  
  expect_s3_class(res, "perm_approx")
  expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
  expect_false(is.null(res$adjust_result))
  
  # lfdr does NOT satisfy monotonicity p_adjust >= p_unadjusted everywhere
  # so we only check that the result is numeric and in [0,1]
  expect_length(res$p_values, length(sim$obs))
})

## ---------------------------------------------------------------------------
## 4. adapt_BH with different true_null_method variants
## ---------------------------------------------------------------------------

test_that("adapt_BH works with all true_null_method options", {
  sim <- simulate_multi(m = 100, n = 100, n_signal = 5)
  
  true_null_methods <- c("convest", "mean", "hist", "farco", "lfdr")
  
  for (mth in true_null_methods) {
    ctrl <- make_adjust_ctrl(true_null_method = mth)
    
    res <- perm_approx(
      obs_stats     = sim$obs,
      perm_stats    = sim$perm,
      method        = "empirical",
      adjust_method = "adapt_BH",
      adjust_ctrl   = ctrl,
      verbose       = FALSE
    )
    
    expect_s3_class(res, "perm_approx")
    expect_false(is.null(res$adjust_result))
    
    # adapt_BH must produce valid adjusted p-values in [0,1]
    expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
    
    # Structure of adjust_result
    expect_true(all(c("p_adjusted", "p_true_null") %in% names(res$adjust_result)))
  }
})

## ---------------------------------------------------------------------------
## 5. Explicit p_true_null overrides estimation
## ---------------------------------------------------------------------------

test_that("p_true_null overrides estimation in adapt_BH", {
  sim <- simulate_multi()
  
  ctrl <- make_adjust_ctrl(true_null_method = "lfdr", p_true_null = 0.6)
  
  res <- perm_approx(
    obs_stats     = sim$obs,
    perm_stats    = sim$perm,
    method        = "empirical",
    adjust_method = "adapt_BH",
    adjust_ctrl   = ctrl,
    verbose       = FALSE
  )
  
  expect_s3_class(res, "perm_approx")
  expect_false(is.null(res$adjust_result))
  
  # p_true_null should be carried through
  expect_equal(res$adjust_result$p_true_null, 0.6, tolerance = 1e-8)
})

## ---------------------------------------------------------------------------
## 6. Empirical + GPD mix under multiple testing correction
##    (tests that adjustment works even when some p-values are GPD-approximated)
## ---------------------------------------------------------------------------

test_that("Adjustment works when mixture of empirical + GPD p-values is present", {
  sim <- simulate_multi(
    m = 20,
    n = 40,
    B = 400,
    n_signal = 4,
    effect = 2.0
  )
  
  gpd_ctrl <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_obs",
    sample_size   = 40,
    exceed0       = 0.2,
    exceed_min    = 5,
    thresh_method = "rob_ftr",
    gof_test      = "none"
  )
  
  res <- perm_approx(
    obs_stats     = sim$obs,
    perm_stats    = sim$perm,
    method        = "gpd",
    gpd_ctrl      = gpd_ctrl,
    approx_thresh = 0.2,
    adjust_method = "BH",
    verbose       = FALSE
  )
  
  expect_s3_class(res, "perm_approx")
  
  # Should have some GPD successes
  expect_true(any(res$method_used == "gpd"))
  
  # BH adjustment must produce valid adjusted p-values
  expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
  
  expect_false(is.null(res$adjust_result))
})

## ---------------------------------------------------------------------------
## 7. Adjustment controls are included in output$control
## ---------------------------------------------------------------------------

test_that("adjust_ctrl appears correctly in output$control", {
  sim <- simulate_multi()
  
  ctrl <- make_adjust_ctrl(true_null_method = "mean")
  
  res <- perm_approx(
    obs_stats     = sim$obs,
    perm_stats    = sim$perm,
    method        = "empirical",
    adjust_method = "BH",
    adjust_ctrl   = ctrl,
    verbose       = FALSE
  )
  
  expect_s3_class(res, "perm_approx")
  
  expect_true("adjust" %in% names(res$control))
  expect_identical(res$control$adjust, ctrl)
})

## ---------------------------------------------------------------------------
## 8. Invalid adjust_ctrl arguments throw expected errors
## ---------------------------------------------------------------------------

test_that("make_adjust_ctrl validation works", {
  # invalid p_true_null
  expect_error(make_adjust_ctrl(p_true_null = -0.1))
  expect_error(make_adjust_ctrl(p_true_null = 1.5))
  expect_error(make_adjust_ctrl(p_true_null = "abc"))
})
