# tests/testthat/test-perm_approx_basic.R

test_that("perm_approx returns valid structure for empirical method", {
  set.seed(42)
  
  n_per_group <- 20
  m_tests     <- 5
  
  group <- rep(c(0, 1), each = n_per_group)
  
  X <- matrix(
    rnorm(2 * n_per_group * m_tests),
    ncol = m_tests
  )
  
  # introduce small effect in first test
  X[group == 1, 1] <- X[group == 1, 1] + 0.5
  
  obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
    colMeans(X[group == 0, , drop = FALSE])
  
  B <- 1000
  perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)
  for (b in seq_len(B)) {
    grp_perm <- sample(group)
    perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
      colMeans(X[grp_perm == 0, , drop = FALSE])
  }
  
  res_emp <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "empirical",
    verbose    = FALSE
  )
  
  # Class and basic structure
  expect_s3_class(res_emp, "perm_approx")
  expect_true(all(c(
    "p_values", "p_unadjusted", "p_empirical",
    "fit_method", "fit_result", "method_used",
    "adjust_result", "control"
  ) %in% names(res_emp)))
  
  # Lengths
  expect_length(res_emp$p_values, m_tests)
  expect_length(res_emp$p_unadjusted, m_tests)
  expect_length(res_emp$p_empirical, m_tests)
  expect_length(res_emp$method_used, m_tests)
  
  # Empirical mode: no parametric fit
  expect_equal(res_emp$fit_method, "empirical")
  expect_null(res_emp$fit_result)
  
  # BH is default: p_values != p_unadjusted in general
  expect_true(any(res_emp$p_values != res_emp$p_unadjusted))
  
  # p-values should be within [0, 1]
  expect_true(all(res_emp$p_unadjusted >= 0 & res_emp$p_unadjusted <= 1))
  expect_true(all(res_emp$p_values      >= 0 & res_emp$p_values      <= 1))
})
