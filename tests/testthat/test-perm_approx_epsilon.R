# tests/testthat/test-perm_approx_epsilon.R

set.seed(123)

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

simulate_perm_for_eps <- function(B = 1000, m = 5) {
  # Simple permutation matrix, roughly N(0,1)
  matrix(rnorm(B * m), nrow = B, ncol = m)
}

.allowed_support_constraints <- c("support_at_obs", "support_at_max")

## ---------------------------------------------------------------------------
## 1. eps_factor: basic behavior
## ---------------------------------------------------------------------------

test_that("eps_factor computes tune * |obs_stats| and ignores perm/sample_size", {
  obs <- c(-2, 0, 3.5)
  perm1 <- matrix(rnorm(100 * 3), nrow = 100)
  perm2 <- matrix(rnorm(100 * 3, mean = 5), nrow = 100)  # very different
  tune  <- 0.1
  
  eps1 <- eps_factor(
    obs_stats   = obs,
    perm_stats  = perm1,
    sample_size = NULL,
    tune        = tune
  )
  
  eps2 <- eps_factor(
    obs_stats   = obs,
    perm_stats  = perm2,
    sample_size = 999,   # arbitrary, should be ignored
    tune        = tune
  )
  
  expected <- tune * abs(obs)
  
  expect_equal(eps1, expected)
  expect_equal(eps2, expected)
})

test_that("eps_factor is monotone in tune", {
  obs  <- c(-1, 0.5, 3)
  perm <- matrix(rnorm(100 * 3), nrow = 100)
  
  tune_small <- 0.05
  tune_big   <- 0.2
  
  eps_small <- eps_factor(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = NULL,
    tune        = tune_small
  )
  
  eps_big <- eps_factor(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = NULL,
    tune        = tune_big
  )
  
  expect_true(all(eps_big >= eps_small))
  # For non-zero obs, strict inequality
  non_zero <- obs != 0
  expect_true(all(eps_big[non_zero] > eps_small[non_zero]))
})

test_that("eps_factor input validation works", {
  perm <- matrix(rnorm(100 * 3), nrow = 100)
  
  expect_error(
    eps_factor(obs_stats = "a", perm_stats = perm, sample_size = NULL, tune = 0.1),
    "must be numeric"
  )
  
  expect_error(
    eps_factor(obs_stats = 1:3, perm_stats = perm, sample_size = NULL, tune = NA),
    "tune' must be a single numeric scalar"
  )
})

## ---------------------------------------------------------------------------
## 2. eps_slls: basic structure and properties
## ---------------------------------------------------------------------------

test_that("eps_slls returns positive, finite eps of correct length", {
  B <- 800
  m <- 4
  perm <- simulate_perm_for_eps(B = B, m = m)
  obs  <- rnorm(m, mean = 2)
  n    <- 50
  
  eps <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = 0.25
  )
  
  expect_type(eps, "double")
  expect_length(eps, length(obs))
  expect_true(all(is.finite(eps)))
  expect_true(all(eps > 0))
})

test_that("eps_slls works with cap_base = 'perm' and cap_base = 'tmax'", {
  B <- 1000
  m <- 5
  perm <- simulate_perm_for_eps(B = B, m = m)
  obs  <- rnorm(m, mean = 1.5)
  n    <- 80
  
  eps_perm <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = 0.25,
    cap_base    = "perm"
  )
  
  eps_tmax <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = 0.25,
    cap_base    = "tmax"
  )
  
  expect_length(eps_perm, length(obs))
  expect_length(eps_tmax, length(obs))
  
  expect_true(all(eps_perm > 0 & is.finite(eps_perm)))
  expect_true(all(eps_tmax > 0 & is.finite(eps_tmax)))
  
  # They don't have to be equal; actually we expect at least some difference
  expect_true(any(abs(eps_perm - eps_tmax) > 0))
})

test_that("eps_slls is monotone in tune (larger tune => larger eps)", {
  B <- 800
  m <- 6
  perm <- simulate_perm_for_eps(B = B, m = m)
  obs  <- rnorm(m, mean = 2)
  n    <- 60
  
  eps_small <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = 0.1
  )
  
  eps_big <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = 0.4
  )
  
  expect_true(all(eps_big >= eps_small - 1e-12))
  expect_true(any(eps_big > eps_small + 1e-8))
})


## ---------------------------------------------------------------------------
## 3. .define_eps: behavior under different constraints & eps_fun
## ---------------------------------------------------------------------------

test_that(".define_eps returns zeros for non-support constraints", {
  perm <- simulate_perm_for_eps(B = 500, m = 4)
  obs  <- rnorm(4)
  
  # unconstrained: eps ignored
  eps_unconstr <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = 50,
    constraint  = "unconstrained",
    eps_fun     = eps_slls,
    eps_tune    = 0.25,
    eps_args    = list()
  )
  
  # shape_nonneg: eps ignored
  eps_shape <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = 50,
    constraint  = "shape_nonneg",
    eps_fun     = eps_slls,
    eps_tune    = 0.25,
    eps_args    = list()
  )
  
  expect_equal(eps_unconstr, rep(0, length(obs)))
  expect_equal(eps_shape,    rep(0, length(obs)))
})

test_that("permApprox:::.define_eps errors when eps_fun is NULL for support constraints", {
  perm <- simulate_perm_for_eps(B = 200, m = 3)
  obs  <- rnorm(3)
  
  expect_error(
    permApprox:::.define_eps(
      perm_stats  = perm,
      obs_stats   = obs,
      sample_size = 50,
      constraint  = "support_at_obs",
      eps_fun     = NULL,
      eps_tune    = 0.25,
      eps_args    = list()
    ),
    "eps_fun' must be non-NULL"
  )
})

test_that("permApprox:::.define_eps with eps_factor matches direct eps_factor call", {
  perm <- simulate_perm_for_eps(B = 400, m = 5)
  obs  <- rnorm(5)
  n    <- 40
  tune <- 0.15
  
  eps_direct <- eps_factor(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = tune
  )
  
  eps_define <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = n,
    constraint  = "support_at_obs",
    eps_fun     = eps_factor,
    eps_tune    = tune,
    eps_args    = list()
  )
  
  expect_equal(eps_define, eps_direct)
})

test_that("permApprox:::.define_eps with eps_slls matches direct eps_slls call (up to numerical noise)", {
  perm <- simulate_perm_for_eps(B = 800, m = 4)
  obs  <- rnorm(4, mean = 1.5)
  n    <- 60
  tune <- 0.3
  
  eps_direct <- eps_slls(
    obs_stats   = obs,
    perm_stats  = perm,
    sample_size = n,
    tune        = tune,
    cap_base    = "perm",
    alpha       = 0.5,
    q_ref       = 0.99,
    rho_lift    = 3,
    k_factor    = 1000,
    floor_epsZ  = 1e-6
  )
  
  eps_define <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = n,
    constraint  = "support_at_obs",
    eps_fun     = eps_slls,
    eps_tune    = tune,
    eps_args    = list(
      cap_base   = "perm",
      alpha      = 0.5,
      q_ref      = 0.99,
      rho_lift   = 3,
      k_factor   = 1000,
      floor_epsZ = 1e-6
    )
  )
  
  expect_equal(eps_define, eps_direct, tolerance = 1e-10)
})

test_that("permApprox:::.define_eps respects eps_args (e.g., cap_base change)", {
  perm <- simulate_perm_for_eps(B = 800, m = 4)
  obs  <- rnorm(4, mean = 1.5)
  n    <- 60
  tune <- 0.25
  
  eps_perm <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = n,
    constraint  = "support_at_obs",
    eps_fun     = eps_slls,
    eps_tune    = tune,
    eps_args    = list(cap_base = "perm")
  )
  
  eps_tmax <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = n,
    constraint  = "support_at_obs",
    eps_fun     = eps_slls,
    eps_tune    = tune,
    eps_args    = list(cap_base = "tmax")
  )
  
  expect_length(eps_perm, length(obs))
  expect_length(eps_tmax, length(obs))
  expect_true(all(eps_perm > 0 & is.finite(eps_perm)))
  expect_true(all(eps_tmax > 0 & is.finite(eps_tmax)))
  
  # eps_args has an effect
  expect_true(any(abs(eps_perm - eps_tmax) > 0))
})

test_that("permApprox:::.define_eps validates the return value of eps_fun", {
  perm <- simulate_perm_for_eps(B = 200, m = 3)
  obs  <- rnorm(3)
  
  bad_fun_len <- function(obs_stats, perm_stats, sample_size, tune, ...) {
    # wrong length
    rep(1, length(obs_stats) + 1)
  }
  
  bad_fun_type <- function(obs_stats, perm_stats, sample_size, tune, ...) {
    # wrong type
    as.character(seq_along(obs_stats))
  }
  
  expect_error(
    permApprox:::.define_eps(
      perm_stats  = perm,
      obs_stats   = obs,
      sample_size = 50,
      constraint  = "support_at_obs",
      eps_fun     = bad_fun_len,
      eps_tune    = 0.25,
      eps_args    = list()
    ),
    "epsilon vector returned by 'eps_fun' must be numeric"
  )
  
  expect_error(
    permApprox:::.define_eps(
      perm_stats  = perm,
      obs_stats   = obs,
      sample_size = 50,
      constraint  = "support_at_obs",
      eps_fun     = bad_fun_type,
      eps_tune    = 0.25,
      eps_args    = list()
    ),
    "epsilon vector returned by 'eps_fun' must be numeric"
  )
})

test_that("permApprox:::.define_eps works with a user-defined eps_fun and eps_args", {
  perm <- simulate_perm_for_eps(B = 300, m = 4)
  obs  <- rnorm(4)
  n    <- 50
  tune <- 0.2
  
  # Custom user-defined eps_fun
  # Uses tune, an extra 'offset' and 'mult' from eps_args
  custom_eps_fun <- function(obs_stats, perm_stats, sample_size, tune,
                             offset = 0, mult = 1, ...) {
    # Just to check the interface: depend on |obs_stats| and the arguments
    base <- abs(obs_stats)
    eps  <- mult * (base + offset) * tune
    eps
  }
  
  # Call permApprox:::.define_eps with constraint = "support_at_obs" and custom eps_fun
  eps_custom <- permApprox:::.define_eps(
    perm_stats  = perm,
    obs_stats   = obs,
    sample_size = n,
    constraint  = "support_at_obs",
    eps_fun     = custom_eps_fun,
    eps_tune    = tune,
    eps_args    = list(offset = 0.5, mult = 2)
  )
  
  # Expected result: 2 * (|obs| + 0.5) * tune
  expected <- 2 * (abs(obs) + 0.5) * tune
  
  expect_type(eps_custom, "double")
  expect_length(eps_custom, length(obs))
  expect_true(all(is.finite(eps_custom)))
  expect_true(all(eps_custom > 0))
  
  # Check that eps_args are actually used
  expect_equal(eps_custom, expected, tolerance = 1e-12)
})
