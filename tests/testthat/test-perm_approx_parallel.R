# tests/testthat/test-perm_approx_basic.R

test_that("perm_approx returns valid structure", {
  set.seed(42)
  
  n_per_group <- 50
  m_tests     <- 14
  B <- 1000
  
  # Group labels: 0 = control, 1 = treatment
  group <- rep(c(0, 1), each = n_per_group)
  
  # Data matrix: rows = samples, cols = tests/features
  X <- matrix(
    rnorm(2 * n_per_group * m_tests),
    ncol = m_tests
  )
  
  ## Introduce true effects:
  ## - Tests 1 and 2: mean(treated) > mean(control)
  ## - Tests 3 and 4: mean(treated) < mean(control)
  true_effects <- c(1.5, 0.8, -0.8, -1.5, rep(0, m_tests - 4))
  
  for (j in seq_len(m_tests)) {
    X[group == 1, j] <- X[group == 1, j] + true_effects[j]
  }
  
  # Observed test statistics: mean difference (treated - control)
  obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
    colMeans(X[group == 0, , drop = FALSE])
  
  # Permutation distribution: shuffle group labels and recompute all 10 stats
  perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)
  
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
      colMeans(X[grp_perm == 0, , drop = FALSE])
  }
  
  # --- GPD approximation ------------------------------------------------------
  gpd_ctrl <- make_gpd_ctrl(
    constraint   = "support_at_obs",
    sample_size  = n_per_group
  )
  
  res_gpd <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gpd",
    alternative = "two_sided",
    gpd_ctrl   = gpd_ctrl,
    verbose    = FALSE,
    cores      = 4
  )
  
  expect_s3_class(res_gpd, "perm_approx")
  expect_equal(res_gpd$fit_method, "gpd")
  expect_false(is.null(res_gpd$fit_result))
  
  fit_gpd <- res_gpd$fit_result
  
  # Expected fields from .compute_pvals_gpd()
  expect_true(all(c(
    "p_value", "thresh", "n_exceed", "shape",
    "scale", "epsilon", "gof_p_value",
    "status", "discrete"
  ) %in% names(fit_gpd)))
  
  # Lengths match number of tests
  expect_length(fit_gpd$p_value, m_tests)
  expect_length(fit_gpd$status,  m_tests)
  
  # GPD p-values used for some tests (success)
  expect_true(any(fit_gpd$status == "success"))
  
  # Final p-values in [0,1]
  expect_true(all(res_gpd$p_unadjusted >= 0 & res_gpd$p_unadjusted <= 1, 
                  na.rm = TRUE))
  
  # --- Gamma approximation ----------------------------------------------------
  gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
  
  res_gamma <- perm_approx(
    obs_stats   = obs_stats,
    perm_stats  = perm_mat,
    method      = "gamma",
    gamma_ctrl  = gamma_ctrl,
    verbose     = FALSE,
    cores       = 4
  )
  
  expect_s3_class(res_gamma, "perm_approx")
  expect_equal(res_gamma$fit_method, "gamma")
  expect_false(is.null(res_gamma$fit_result))
  
  fit_gamma <- res_gamma$fit_result
  
  expect_true(all(c(
    "p_value", "shape", "rate",
    "gof_p_value", "status", "discrete",
    "gof_rejected", "method_used"
  ) %in% names(fit_gamma)))
  
  expect_length(fit_gamma$p_value, m_tests)
  expect_length(fit_gamma$status,  m_tests)
  
  # Success for at least some tests
  expect_true(any(fit_gamma$status == "success"))
  
  # Final p-values in [0,1]
  expect_true(all(res_gamma$p_unadjusted >= 0 & res_gamma$p_unadjusted <= 1, 
                  na.rm = TRUE))
})
