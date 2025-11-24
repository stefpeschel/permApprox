# tests/testthat/test-gpd_fit_variations.R

set.seed(42)

## ---------------------------------------------------------------------------
## Helper: simple two-group simulation used across tests
## ---------------------------------------------------------------------------

simulate_two_group <- function(
    n_per_group = 40,
    m_tests     = 8,
    effects     = NULL,
    B           = 500
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

.allowed_gpd_status <- c(
  "not_selected", "success", "discrete",
  "no_threshold", "gof_reject", "fit_failed"
)

## ---------------------------------------------------------------------------
## 1. Basic GPD LME fit under support_at_obs constraint
## ---------------------------------------------------------------------------

test_that("GPD LME with support_at_obs returns valid structure and sane p-values", {
  n_per_group <- 50
  sim <- simulate_two_group(
    n_per_group = n_per_group,
    m_tests     = 10,
    effects     = c(1.5, 1, 0.8, rep(0, 7)),
    B           = 800
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  gpd_ctrl <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_obs",
    sample_size   = n_per_group,            # required for eps_slls
    thresh_method = "rob_ftr",
    exceed0       = 0.25,
    exceed_min    = 10,
    gof_test      = "none"
  )
  
  res_gpd <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.2,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_ctrl,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  expect_s3_class(res_gpd, "perm_approx")
  expect_identical(res_gpd$fit_method, "gpd")
  
  fit <- res_gpd$fit_result
  
  # Expected fields from .compute_pvals_gpd()
  expect_true(all(c(
    "p_value", "thresh", "n_exceed", "shape",
    "scale", "epsilon", "gof_p_value",
    "status", "discrete"
  ) %in% names(fit)))
  
  m <- length(obs)
  
  # Lengths
  expect_length(fit$p_value,   m)
  expect_length(fit$thresh,    m)
  expect_length(fit$n_exceed,  m)
  expect_length(fit$shape,     m)
  expect_length(fit$scale,     m)
  expect_length(fit$epsilon,   m)
  expect_length(fit$status,    m)
  expect_length(fit$discrete,  m)
  
  # Status in allowed set
  expect_true(all(as.character(fit$status) %in% .allowed_gpd_status))
  
  # p-values in [0,1]
  expect_true(all(res_gpd$p_unadjusted >= 0 & res_gpd$p_unadjusted <= 1, 
                  na.rm = TRUE))
  
  # For successful GPD fits, check some sanity
  idx_success <- which(fit$status == "success")
  
  expect_true(all(is.finite(fit$shape[idx_success])))
  expect_true(all(is.finite(fit$scale[idx_success]) & fit$scale[idx_success] > 0))
  
  # support constraint: epsilon should be non-negative and finite
  expect_true(all(fit$epsilon[idx_success] >= 0))
  expect_true(all(is.finite(fit$epsilon[idx_success])))
  
  # ensure some GPD p-values differ from empirical
  expect_true(any(res_gpd$p_unadjusted[idx_success] != res_gpd$p_empirical[idx_success]))
  
  # No GPD-based p-value should be exactly 0 under support constraint
  expect_true(all(res_gpd$p_unadjusted[idx_success] > 0))
})

## ---------------------------------------------------------------------------
## 2. All GPD fit_method options (unconstrained) run and give valid p-values
## ---------------------------------------------------------------------------

test_that("All GPD fit_method options run (unconstrained) and produce sensible outputs", {
  sim <- simulate_two_group(
    n_per_group = 40,
    m_tests     = 6,
    effects     = c(1.5, 1, rep(0, 4)),  # a few signals
    B           = 600
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  fit_methods <- c("lme", "mle1d", "mle2d", "mom", "nls2", "wnllsm", "zse")
  
  for (fm in fit_methods) {
    gpd_ctrl <- make_gpd_ctrl(
      fit_method    = fm,
      constraint    = "unconstrained",
      thresh_method = "ftr",
      exceed0       = 0.2,
      exceed_min    = 5,
      gof_test      = "none"
    )
    
    res <- perm_approx(
      obs_stats      = obs,
      perm_stats     = perm,
      method         = "gpd",
      approx_thresh  = 0.5,          # encourage GPD use
      alternative    = "two_sided",
      gpd_ctrl       = gpd_ctrl,
      adjust_method  = "none",
      verbose        = FALSE
    )
    
    fit <- res$fit_result
    
    expect_s3_class(res, "perm_approx")
    expect_identical(res$fit_method, "gpd")
    
    # p-values in [0,1]
    expect_true(all(res$p_unadjusted >= 0 & res$p_unadjusted <= 1, na.rm = TRUE))
    
    # structure
    expect_true(all(c(
      "p_value", "thresh", "n_exceed", "shape",
      "scale", "epsilon", "gof_p_value",
      "status", "discrete"
    ) %in% names(fit)))
    
    expect_true(all(as.character(fit$status) %in% .allowed_gpd_status))
    
    # At least one successful fit for most methods (if not, we just skip)
    if (any(fit$status == "success")) {
      idx_success <- which(fit$status == "success")
      expect_true(all(is.finite(fit$shape[idx_success])))
      expect_true(all(is.finite(fit$scale[idx_success]) & fit$scale[idx_success] > 0))
    }
  }
})

## ---------------------------------------------------------------------------
## 3. shape_nonneg constraint enforces non-negative shape (mle1d)
## ---------------------------------------------------------------------------

test_that("shape_nonneg constraint with mle1d yields non-negative shape", {
  sim <- simulate_two_group(
    n_per_group = 50,
    m_tests     = 6,
    effects     = c(1.5, 2, rep(0, 4)),
    B           = 700
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  gpd_ctrl <- make_gpd_ctrl(
    fit_method    = "mle1d",
    constraint    = "shape_nonneg",
    thresh_method = "ftr",
    exceed0       = 0.3,
    exceed_min    = 5,
    gof_test      = "none"
  )
  
  res <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.4,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_ctrl,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  fit <- res$fit_result
  
  expect_s3_class(res, "perm_approx")
  expect_identical(res$fit_method, "gpd")
  
  idx_success <- which(fit$status == "success")
  if (length(idx_success) > 0L) {
    # allow tiny numerical negatives, but not substantial ones
    expect_true(all(fit$shape[idx_success] > -1e-6))
  }
  
  # No zero p-values
  expect_true(all(res$p_unadjusted[idx_success] > 0))
})

## ---------------------------------------------------------------------------
## 4. support_at_obs vs support_at_max constraints
## ---------------------------------------------------------------------------

test_that("support_at_obs and support_at_max both yield non-zero GPD p-values", {
  n_per_group <- 60
  sim <- simulate_two_group(
    n_per_group = n_per_group,
    m_tests     = 12,
    effects     = c(1.8, 1.2, 1, rep(0, 9)),
    B           = 800
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  # Per-test support constraint
  gpd_obs <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_obs",
    sample_size   = n_per_group,
    thresh_method = "rob_ftr",
    exceed0       = 0.2,
    exceed_min    = 10,
    gof_test      = "none",
    zero_guard    = TRUE
  )
  
  # Global support constraint at max observed statistic
  gpd_max <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_max",
    sample_size   = n_per_group,
    thresh_method = "rob_ftr",
    exceed0       = 0.2,
    exceed_min    = 10,
    gof_test      = "none",
    zero_guard    = TRUE
  )
  
  res_obs <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.2,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_obs,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  res_max <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.2,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_max,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  fit_obs <- res_obs$fit_result
  fit_max <- res_max$fit_result
  
  idx_succ_obs <- which(fit_obs$status == "success")
  idx_succ_max <- which(fit_max$status == "success")
  
  # For successful GPD fits, no zero p-values
  expect_true(all(res_obs$p_unadjusted[idx_succ_obs] > 0))
  expect_true(all(res_max$p_unadjusted[idx_succ_max] > 0))
  
  # epsilon is defined and finite for successful fits
  expect_true(all(fit_obs$epsilon[idx_succ_obs] >= 0))
  expect_true(all(is.finite(fit_obs$epsilon[idx_succ_obs])))
  expect_true(all(fit_max$epsilon[idx_succ_max] >= 0))
  expect_true(all(is.finite(fit_max$epsilon[idx_succ_max])))
})

## ---------------------------------------------------------------------------
## 5. Zero-guard on vs off both run and give valid results
## ---------------------------------------------------------------------------

test_that("zero_guard on/off both run and yield valid GPD outputs", {
  n_per_group <- 5000
  sim <- simulate_two_group(
    n_per_group = n_per_group,
    m_tests     = 10,
    effects     = c(50, 1.5, 1, rep(0, 7)),  # strong signals
    B           = 1000
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  gpd_zg_on <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_obs",
    sample_size   = n_per_group,
    thresh_method = "rob_ftr",
    exceed0       = 0.2,
    exceed_min    = 10,
    gof_test      = "none",
    zero_guard    = TRUE
  )
  
  gpd_zg_off <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "support_at_obs",
    sample_size   = n_per_group,
    thresh_method = "rob_ftr",
    exceed0       = 0.2,
    exceed_min    = 10,
    gof_test      = "none",
    zero_guard    = FALSE
  )
  
  res_on <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.1,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_zg_on,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  res_off <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.1,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_zg_off,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  # Only non-zero p-values with zero guard
  expect_true(all(res_on$p_unadjusted > 0))
  
  # Structural checks
  for (res in list(res_on, res_off)) {
    expect_s3_class(res, "perm_approx")
    expect_identical(res$fit_method, "gpd")
    expect_true(all(res$p_unadjusted >= 0 & res$p_unadjusted <= 1, na.rm = TRUE))
    
    fit <- res$fit_result
    expect_true(all(as.character(fit$status) %in% .allowed_gpd_status))
    
    idx_succ <- which(fit$status == "success")
    expect_true(all(fit$epsilon[idx_succ] >= 0))
    expect_true(all(is.finite(fit$epsilon[idx_succ])))
    expect_true(all(is.finite(fit$shape[idx_succ])))
    expect_true(all(is.finite(fit$scale[idx_succ]) & fit$scale[idx_succ] > 0))
    
  }
})

## ---------------------------------------------------------------------------
## 6. Threshold methods: 'fix' vs default 'rob_ftr'
## ---------------------------------------------------------------------------

test_that("GPD works with thresh_method = 'fix' as well as 'rob_ftr'", {
  sim <- simulate_two_group(
    n_per_group = 50,
    m_tests     = 8,
    effects     = c(1.5, 1, rep(0, 6)),
    B           = 600
  )
  obs  <- sim$obs
  perm <- sim$perm_mat
  
  # low fixed threshold to guarantee many exceedances (even if not ideal in practice)
  fixed_thresh <- min(perm) - 1
  
  gpd_fix <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "unconstrained",
    thresh_method = "fix",
    thresh0       = fixed_thresh,
    exceed0       = NULL,
    exceed_min    = 10,
    gof_test      = "none"
  )
  
  gpd_rob <- make_gpd_ctrl(
    fit_method    = "lme",
    constraint    = "unconstrained",
    thresh_method = "rob_ftr",
    exceed0       = 0.2,
    exceed_min    = 10,
    gof_test      = "none"
  )
  
  res_fix <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.2,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_fix,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  res_rob <- perm_approx(
    obs_stats      = obs,
    perm_stats     = perm,
    method         = "gpd",
    approx_thresh  = 0.2,
    alternative    = "two_sided",
    gpd_ctrl       = gpd_rob,
    adjust_method  = "none",
    verbose        = FALSE
  )
  
  for (res in list(res_fix, res_rob)) {
    fit <- res$fit_result
    
    expect_s3_class(res, "perm_approx")
    expect_identical(res$fit_method, "gpd")
    expect_true(all(res$p_unadjusted >= 0 & res$p_unadjusted <= 1, na.rm = TRUE))
    expect_true(all(as.character(fit$status) %in% .allowed_gpd_status))
    
    idx_succ <- which(fit$status == "success")
    if (length(idx_succ) > 0L) {
      expect_true(all(fit$n_exceed[idx_succ] >= gpd_fix$exceed_min))
    }
  }
})
