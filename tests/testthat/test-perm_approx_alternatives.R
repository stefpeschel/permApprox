# tests/testthat/test-perm_approx_alternatives.R

test_that("perm_approx works all alternatives (two_sided, less, greater)", {
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
  
  ##############################################################################
  ### Empirical p-values
  ##############################################################################
  
  # Compute empirical p-values
  res_emp_greater <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "empirical",
    alternative = "greater",
    verbose    = FALSE
  )
  
  res_emp_less <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "empirical",
    alternative = "less",
    verbose    = FALSE
  )
  
  res_emp_2sided <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "empirical",
    alternative = "two_sided",
    verbose    = FALSE
  )
  
  # P-values of significant tests correct?
  expect_true(all(res_emp_greater$p_unadjusted[1:2] == 1/(B+1)))
  expect_true(all(res_emp_less$p_unadjusted[3:4] == 1/(B+1)))
  expect_true(all(res_emp_2sided$p_unadjusted[1:4] == 1/(B+1)))
  
  # p-values should be within [0, 1]
  expect_true(all(res_emp_greater$p_unadjusted >= 0 & 
                    res_emp_greater$p_unadjusted <= 1))
  expect_true(all(res_emp_less$p_unadjusted >= 0 & 
                    res_emp_less$p_unadjusted <= 1))
  expect_true(all(res_emp_2sided$p_unadjusted >= 0 & 
                    res_emp_2sided$p_unadjusted <= 1))
  
  ##############################################################################
  ### GPD approximation
  ##############################################################################
  
  # --- Two-sided --------------------------------------------------------------
  
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
    verbose    = FALSE
  )
  
  # Final p-values in [0,1]
  expect_true(all(res_gpd$p_unadjusted >= 0 & res_gpd$p_unadjusted <= 1, 
                  na.rm = TRUE))
  
  expect_true(all(res_gpd$p_unadjusted[1:4] < 0.01))
  expect_true(all(res_gpd$p_unadjusted[5:14] > 0.01))
  expect_true(res_gpd$p_unadjusted[1] < res_gpd$p_unadjusted[2])
  expect_true(res_gpd$p_unadjusted[4] < res_gpd$p_unadjusted[3])
  
  ## --- Right-sided -----------------------------------------------------------
  res_gpd_greater <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gpd",
    alternative = "greater",
    gpd_ctrl   = gpd_ctrl,
    verbose    = FALSE
  )
  
  expect_true(all(res_gpd_greater$p_unadjusted[1:2] < 0.01))
  expect_true(all(res_gpd_greater$p_unadjusted[3:4] == 1))
  expect_true(all(res_gpd_greater$p_unadjusted[3:14] > 0.01))
  expect_true(res_gpd$p_unadjusted[1] < res_gpd$p_unadjusted[2])
  
  ## --- Left-sided ------------------------------------------------------------
  res_gpd_less <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gpd",
    alternative = "less",
    gpd_ctrl   = gpd_ctrl,
    verbose    = FALSE
  )
  
  expect_true(all(res_gpd_less$p_unadjusted[3:4] < 0.01))
  expect_true(all(res_gpd_less$p_unadjusted[1:2] == 1))
  expect_true(all(res_gpd_less$p_unadjusted[5:14] > 0.01))
  expect_true(res_gpd$p_unadjusted[4] < res_gpd$p_unadjusted[3])
  
  ##############################################################################
  ### Gamma approximation
  ##############################################################################
  
  # --- Two-sided --------------------------------------------------------------
  
  gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
  
  res_gamma <- perm_approx(
    obs_stats   = obs_stats,
    perm_stats  = perm_mat,
    method      = "gamma",
    gamma_ctrl  = gamma_ctrl,
    verbose     = FALSE
  )
  
  # Final p-values in [0,1]
  expect_true(all(res_gamma$p_unadjusted >= 0 & res_gamma$p_unadjusted <= 1, 
                  na.rm = TRUE))
  
  expect_true(all(res_gamma$p_unadjusted[1:4] < 0.01))
  expect_true(all(res_gamma$p_unadjusted[5:14] > 0.01))
  expect_true(res_gamma$p_unadjusted[1] < res_gamma$p_unadjusted[2])
  expect_true(res_gamma$p_unadjusted[4] < res_gamma$p_unadjusted[3])
  
  ## --- Right-sided -----------------------------------------------------------
  res_gamma_greater <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gamma",
    alternative = "greater",
    gamma_ctrl  = gamma_ctrl,
    verbose    = FALSE
  )
  
  expect_true(all(res_gamma_greater$p_unadjusted[1:2] < 0.01))
  expect_true(all(res_gamma_greater$p_unadjusted[3:4] == 1))
  expect_true(all(res_gamma_greater$p_unadjusted[3:14] > 0.01))
  expect_true(res_gamma$p_unadjusted[1] < res_gamma$p_unadjusted[2])
  
  ## --- Left-sided ------------------------------------------------------------
  res_gamma_less <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gamma",
    alternative = "less",
    gamma_ctrl  = gamma_ctrl,
    verbose    = FALSE
  )
  
  expect_true(all(res_gamma_less$p_unadjusted[3:4] < 0.01))
  expect_true(all(res_gamma_less$p_unadjusted[1:2] == 1))
  expect_true(all(res_gamma_less$p_unadjusted[5:14] > 0.01))
  expect_true(res_gamma$p_unadjusted[4] < res_gamma$p_unadjusted[3])
})
