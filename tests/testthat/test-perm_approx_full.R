# tests/testthat/test-perm_approx_full.R

set.seed(42)

## Helper to simulate a simple two-group setup -------------------------------

simulate_two_group <- function(
    n_per_group = 30,
    m_tests     = 8,
    effects     = NULL,
    B           = 300
) {
  if (is.null(effects)) {
    effects <- rep(0, m_tests)
  }
  stopifnot(length(effects) == m_tests)
  
  group <- rep(c(0, 1), each = n_per_group)
  
  X <- matrix(
    rnorm(2 * n_per_group * m_tests),
    ncol = m_tests
  )
  
  for (j in seq_len(m_tests)) {
    X[group == 1, j] <- X[group == 1, j] + effects[j]
  }
  
  obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
    colMeans(X[group == 0, , drop = FALSE])
  
  perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
      colMeans(X[grp_perm == 0, , drop = FALSE])
  }
  
  list(
    group    = group,
    X        = X,
    obs      = obs_stats,
    perm_mat = perm_mat
  )
}

## ---------------------------------------------------------------------------
## 1. Argument validation and basic invariants
## ---------------------------------------------------------------------------

test_that("perm_approx does basic argument validation", {
  gpd_ctrl <- make_gpd_ctrl(constraint = "unconstrained")
  
  sim <- simulate_two_group(m_tests = 3)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  # obs_stats length vs perm_stats ncol
  expect_error(
    perm_approx(obs_stats = obs[1:2], perm_stats = perm, gpd_ctrl = gpd_ctrl),
    "Length of 'obs_stats' must match number of columns"
  )
  
  # invalid alternative
  expect_error(
    perm_approx(
      obs_stats  = obs,
      perm_stats = perm,
      gpd_ctrl = gpd_ctrl,
      alternative = "two.sided"  # wrong spelling
    ),
    "'arg' should be one of"
  )
  
  # invalid method
  expect_error(
    perm_approx(
      obs_stats  = obs,
      perm_stats = perm,
      gpd_ctrl = gpd_ctrl,
      method     = "gpds"
    ),
    "'arg' should be one of"
  )
  
  # invalid power
  expect_error(
    perm_approx(
      obs_stats  = obs,
      perm_stats = perm,
      gpd_ctrl = gpd_ctrl,
      power      = 0
    ),
    "power > 0"
  )
  
  # invalid cores
  expect_error(
    perm_approx(
      obs_stats  = obs,
      perm_stats = perm,
      gpd_ctrl = gpd_ctrl,
      cores      = 0
    ),
    "must be a single integer >= 1"
  )
  
  # invalid parallel_min
  expect_error(
    perm_approx(
      obs_stats    = obs,
      perm_stats   = perm,
      gpd_ctrl = gpd_ctrl,
      parallel_min = 0L
    ),
    "must be a single integer >= 1"
  )
  
  # invalid verbose
  expect_error(
    perm_approx(
      obs_stats  = obs,
      perm_stats = perm,
      gpd_ctrl = gpd_ctrl,
      verbose    = NA
    ),
    "'verbose' must be a single logical"
  )
  
  # invalid adjust_ctrl
  expect_error(
    perm_approx(
      obs_stats   = obs,
      perm_stats  = perm,
      gpd_ctrl = gpd_ctrl,
      adjust_ctrl = list()
    ),
    "must be created with make_adjust_ctrl"
  )
})

## ---------------------------------------------------------------------------
## 2. Single-test usage with vector perm_stats
## ---------------------------------------------------------------------------

test_that("perm_approx works with single test and vector perm_stats", {
  set.seed(1)
  B <- 200
  n_per_group <- 20
  
  # Null situation: one test, no effect
  group <- rep(c(0, 1), each = n_per_group)
  x     <- rnorm(n_per_group * 2)
  
  obs <- mean(x[group == 1]) - mean(x[group == 0])
  
  perm_vec <- numeric(B)
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm_vec[b] <- mean(x[grp_perm == 1]) - mean(x[grp_perm == 0])
  }
  
  # Use vector for perm_stats + scalar obs_stats
  res <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm_vec,
    method         = "empirical",
    alternative    = "greater",
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  expect_s3_class(res, "perm_approx")
  expect_length(res$p_values, 1L)
  expect_length(res$p_unadjusted, 1L)
  expect_length(res$p_empirical, 1L)
  expect_length(res$n_perm_exceed, 1L)
  
  # Check manual empirical p-value
  r_manual <- sum(perm_vec >= obs)
  p_manual <- (r_manual + 1) / (B + 1)
  
  expect_equal(as.numeric(res$p_unadjusted), p_manual)
  expect_equal(as.numeric(res$p_empirical), p_manual)
  expect_identical(res$p_values, res$p_unadjusted) # no adjustment
})

## ---------------------------------------------------------------------------
## 3. Centering: numeric vs "mean" / "median"
## ---------------------------------------------------------------------------

test_that("null_center='mean' equals numeric centering at column means", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 4, n_per_group = n_per_group)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  col_means <- colMeans(perm)
  
  res_mean_keyword <- perm_approx(
    obs_stats   = obs,
    perm_stats  = perm,
    method      = "empirical",
    null_center = "mean",
    adjust_method = "none",
    verbose     = FALSE
  )
  
  res_mean_numeric <- perm_approx(
    obs_stats   = obs,
    perm_stats  = perm,
    method      = "empirical",
    null_center = col_means,
    adjust_method = "none",
    verbose     = FALSE
  )
  
  expect_equal(res_mean_keyword$p_empirical, res_mean_numeric$p_empirical)
  expect_equal(res_mean_keyword$p_unadjusted, res_mean_numeric$p_unadjusted)
})

test_that("null_center='median' equals numeric centering at column medians", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 4, n_per_group = n_per_group)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  col_medians <- apply(perm, 2L, median)
  
  res_med_keyword <- perm_approx(
    obs_stats   = obs,
    perm_stats  = perm,
    method      = "empirical",
    null_center = "median",
    adjust_method = "none",
    verbose     = FALSE
  )
  
  res_med_numeric <- perm_approx(
    obs_stats   = obs,
    perm_stats  = perm,
    method      = "empirical",
    null_center = col_medians,
    adjust_method = "none",
    verbose     = FALSE
  )
  
  expect_equal(res_med_keyword$p_empirical, res_med_numeric$p_empirical)
  expect_equal(res_med_keyword$p_unadjusted, res_med_numeric$p_unadjusted)
})

## ---------------------------------------------------------------------------
## 4. Power transformation consistency
## ---------------------------------------------------------------------------

test_that("power transformation matches manual transformation", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 3, B = 200, n_per_group = n_per_group)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  power <- 2
  
  # Run perm_approx with power=2
  res_power <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "empirical",
    power          = power,
    alternative    = "greater",
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  # Manual: transform first, then compute empirical p-values
  obs_tr   <- obs^power
  perm_tr  <- perm^power
  
  r_manual <- vapply(
    seq_along(obs_tr),
    function(j) sum(perm_tr[, j] >= obs_tr[j]),
    integer(1)
  )
  p_manual <- (r_manual + 1) / (nrow(perm_tr) + 1)
  
  expect_equal(res_power$p_empirical, p_manual)
  expect_equal(res_power$p_unadjusted, p_manual)
})

## ---------------------------------------------------------------------------
## 5. No tests selected for parametric approximation (idx_fit empty)
## ---------------------------------------------------------------------------

test_that("parametric methods fall back to empirical when no test is selected", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 5, B = 300, n_per_group = n_per_group)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  # Baseline empirical result
  res_emp <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "empirical",
    approx_thresh  = 0.1,
    alternative    = "two_sided",
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  # approx_thresh = 0 => p_empirical < approx_thresh can never hold
  # => idx_fit is empty
  gpd_ctrl <- make_gpd_ctrl(sample_size = n_per_group)
  gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
  
  res_gpd <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_ctrl,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  res_gamma <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gamma",
    approx_thresh  = 0,
    alternative    = "two_sided",
    gamma_ctrl     = gamma_ctrl,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  # In both cases no test should be replaced by parametric approximation
  expect_equal(res_gpd$p_unadjusted,   res_emp$p_unadjusted)
  expect_equal(res_gamma$p_unadjusted, res_emp$p_unadjusted)
  
  expect_true(all(res_gpd$method_used   == "empirical"))
  expect_true(all(res_gamma$method_used == "empirical"))
  
  expect_identical(res_gpd$fit_method, "gpd")
  expect_identical(res_gamma$fit_method, "gamma")
  
  # fit_result should exist, but contain no successful fits
  expect_false(is.null(res_gpd$fit_result))
  expect_false(is.null(res_gamma$fit_result))
  
  expect_true(all(res_gpd$fit_result$status %in% 
                    c("not_selected", "no_threshold", "fit_failed", "discrete", "gof_reject")))
  expect_true(all(res_gamma$fit_result$status %in% 
                    c("not_selected", "fit_failed", "discrete", "gof_reject")))
})

## ---------------------------------------------------------------------------
## 6. Multiple testing adjustment methods
## ---------------------------------------------------------------------------

test_that("adjust_method = 'none' leaves p_values unchanged", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 6, B = 300, n_per_group = n_per_group)
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  res_none <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    gpd_ctrl       = make_gpd_ctrl(constraint = "support_at_obs", 
                                   sample_size = n_per_group),
    approx_thresh  = 0.2,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  expect_true(all(res_none$p_values == res_none$p_unadjusted))
  expect_null(res_none$adjust_result)
})

test_that("adjust_method from p.adjust (e.g. 'holm') works", {
  sim <- simulate_two_group(m_tests = 10, B = 300,
                            effects = c(1.5, 1, rep(0, 8)))
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  res_holm <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "empirical",
    adjust_method  = "holm",
    verbose        = FALSE
  )
  
  expect_s3_class(res_holm, "perm_approx")
  expect_true(all(res_holm$p_values >= res_holm$p_unadjusted))
  expect_true(all(res_holm$p_values >= 0 & res_holm$p_values <= 1))
})

test_that("adapt_BH adjustment runs and returns valid p-values", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 12, B = 300, n_per_group = n_per_group,
                            effects = c(1.5, 1, 0.8, rep(0, 9)))
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  gpd_ctrl    <- make_gpd_ctrl(constraint = "support_at_obs", 
                               sample_size = n_per_group)
  adjust_ctrl <- make_adjust_ctrl(true_null_method = "lfdr")
  
  res_adapt <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    gpd_ctrl       = gpd_ctrl,
    adjust_method  = "adapt_BH",
    adjust_ctrl    = adjust_ctrl,
    verbose        = FALSE
  )
  
  expect_s3_class(res_adapt, "perm_approx")
  expect_equal(length(res_adapt$p_values), length(obs))
  expect_true(all(res_adapt$p_values >= 0 & res_adapt$p_values <= 1, na.rm = TRUE))
})

## ---------------------------------------------------------------------------
## 7. Alternative handling and transformed sign constraint
## ---------------------------------------------------------------------------

test_that("tests with non-positive transformed obs_stat are not approximated", {
  n_per_group <- 30
  sim <- simulate_two_group(
    m_tests = 4,
    effects = c(-1.5, -1, 0, 0),    # negative effects for first two
    B       = 300,
    n_per_group = n_per_group
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  # 'greater' alternative: negative obs_stats should map to t_obs <= 0
  # and thus be excluded from GPD/Gamma fitting.
  gpd_ctrl <- make_gpd_ctrl(constraint = "support_at_obs", 
                            sample_size = n_per_group)
  
  res_gpd <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    alternative    = "greater",
    approx_thresh  = 0.5,
    gpd_ctrl       = gpd_ctrl,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  # Some tests may still be empirical due to GOF issues, but
  # negative observed statistics (indices 1:2) must not use GPD.
  expect_true(all(res_gpd$method_used[1:2] == "empirical"))
})

## ---------------------------------------------------------------------------
## 8. Print and summary methods do not error
## ---------------------------------------------------------------------------

test_that("print() and summary() for perm_approx run without error", {
  n_per_group <- 30
  sim <- simulate_two_group(m_tests = 5, B = 1000, n_per_group = n_per_group,
                            effects = c(1.5, rep(0, 4)))
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  res <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    gpd_ctrl       = make_gpd_ctrl(constraint = "support_at_obs",
                                   sample_size = n_per_group),
    verbose        = FALSE
  )
  
  expect_silent(capture.output(print(res)))
  expect_silent(capture.output(summary(res)))
})
