# tests/testthat/test-print-summary.R

test_that("print.perm_approx and summary.perm_approx run without error", {
  set.seed(42)
  
  n_per_group <- 20
  m_tests     <- 3
  
  group <- rep(c(0, 1), each = n_per_group)
  
  X <- matrix(
    rnorm(2 * n_per_group * m_tests),
    ncol = m_tests
  )
  
  # Add small effect in first test
  X[group == 1, 1] <- X[group == 1, 1] + 0.6
  
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
  
  res <- perm_approx(
    obs_stats  = obs_stats,
    perm_stats = perm_mat,
    method     = "gpd",
    gpd_ctrl   = gpd_ctrl,
    verbose    = FALSE
  )
  
  # print() should not error and should return object invisibly
  out_print <- testthat::capture_output(
    invisible_res <- print(res)
  )
  expect_s3_class(invisible_res, "perm_approx")
  expect_true(nchar(out_print) > 0)
  
  # summary() likewise
  out_summary <- testthat::capture_output(
    invisible_sum <- summary(res)
  )
  expect_s3_class(invisible_sum, "perm_approx")
  expect_true(nchar(out_summary) > 0)
})
