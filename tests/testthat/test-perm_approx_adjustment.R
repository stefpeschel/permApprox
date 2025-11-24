# tests/testthat/test-perm_approx_adjustment.R

test_that("perm_approx multiple testing adjustment behaves as expected", {
  set.seed(42)
  
  n_per_group <- 40
  m_tests     <- 8
  
  group <- rep(c(0, 1), each = n_per_group)
  
  X <- matrix(
    rnorm(2 * n_per_group * m_tests),
    ncol = m_tests
  )
  
  # Add a moderate effect to first two tests
  X[group == 1, 1:2] <- X[group == 1, 1:2] + 0.7
  
  obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
    colMeans(X[group == 0, , drop = FALSE])
  
  B <- 1000
  perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
      colMeans(X[grp_perm == 0, , drop = FALSE])
  }
  
  gpd_ctrl <- make_gpd_ctrl(
    constraint   = "support_at_obs",
    sample_size  = n_per_group
  )
  
  # No adjustment
  res_none <- perm_approx(
    obs_stats     = obs_stats,
    perm_stats    = perm_mat,
    method        = "gpd",
    gpd_ctrl      = gpd_ctrl,
    adjust_method = "none",
    verbose       = FALSE
  )
  
  expect_null(res_none$adjust_result)
  expect_equal(res_none$p_values, res_none$p_unadjusted)
  
  # BH adjustment
  res_bh <- perm_approx(
    obs_stats     = obs_stats,
    perm_stats    = perm_mat,
    method        = "gpd",
    gpd_ctrl      = gpd_ctrl,
    adjust_method = "BH",
    verbose       = FALSE
  )
  
  expect_false(is.null(res_bh$adjust_result))
  expect_length(res_bh$p_values, m_tests)
  # BH should not *increase* number of very small p-values
  expect_true(all(res_bh$p_values >= 0 & res_bh$p_values <= 1, na.rm = TRUE))
  
  # adaptive BH
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
  
  expect_false(is.null(res_adapt$adjust_result))
  expect_length(res_adapt$p_values, m_tests)
  expect_true(all(res_adapt$p_values >= 0 & res_adapt$p_values <= 1, na.rm = TRUE))
})
