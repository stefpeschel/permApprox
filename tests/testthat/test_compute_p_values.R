library(testthat)
library(permAprox)

test_that("compute_p_values runs without errors for all argument combinations", {
  set.seed(123)
  # Single and multiple test setups
  obs_vec <- 2.0
  perm_vec <- rnorm(1000)
  obs_mat <- c(2.0, 3.0, 4.0, 5.0)
  perm_mat <- matrix(rnorm(4 * 1000), nrow = 4)

  # Argument options
  methods <- c("empirical", "gamma", "gpd")
  alternatives <- c("greater", "less", "two_sided")
  null_centers <- list(0, "mean", "median")
  adjust_methods <- c("none", "BH")

  # Control objects
  gpd_ctrl <- make_gpd_ctrl()
  gamma_ctrl <- make_gamma_ctrl()
  adjust_ctrl <- make_adjust_ctrl()

  # Test all combinations
  for (method in methods) {
    for (alternative in alternatives) {
      for (null_center in null_centers) {
        for (adjust_method in adjust_methods) {
          # Multiple tests (matrix input)
          expect_error(
            compute_p_values(
              obs_stats = obs_mat,
              perm_stats = perm_mat,
              method = method,
              fit_thresh = 0.2,
              alternative = alternative,
              null_center = null_center,
              adjust_method = adjust_method,
              gpd_ctrl = gpd_ctrl,
              gamma_ctrl = gamma_ctrl,
              adjust_ctrl = adjust_ctrl
            ),
            NA
          )
          # Single test (vector input)
          expect_error(
            compute_p_values(
              obs_stats = obs_vec,
              perm_stats = perm_vec,
              method = method,
              fit_thresh = 0.2,
              alternative = alternative,
              null_center = null_center,
              adjust_method = adjust_method,
              gpd_ctrl = gpd_ctrl,
              gamma_ctrl = gamma_ctrl,
              adjust_ctrl = adjust_ctrl
            ),
            NA
          )
        }
      }
    }
  }
})


# Test compute_p_values with custom control parameters

test_that("compute_p_values accepts custom control parameters without errors", {
  set.seed(123)
  obs_vec <- 2.0
  perm_vec <- rnorm(1000)
  obs_mat <- c(2.0, 3.0, 4.0, 5.0)
  perm_mat <- matrix(rnorm(4 * 1000), nrow = 4)

  # Custom GPD control (non-default fit method and include_obs)
  gpd_ctrl_custom <- make_gpd_ctrl(fit_method = "mom",
                                   include_obs = TRUE)

  # Custom Gamma control (include observed, no GOF test, different alpha)
  gamma_ctrl_custom <- make_gamma_ctrl(include_obs = TRUE,
                                       gof_test = "none",
                                       gof_alpha = 0.1)

  # Custom adjust control (empirical true null, p_true_null set)
  adjust_ctrl_custom <- make_adjust_ctrl(true_null_method = "farco",
                                         p_true_null = 0.2,
                                         seq_length = 50,
                                         cores = 2,
                                         verbose = TRUE)

  # A representative parameter set
  expect_error(
    compute_p_values(
      obs_stats = obs_vec,
      perm_stats = perm_vec,
      method = "gpd",
      fit_thresh = 0.3,
      alternative = "greater",
      null_center = "median",
      adjust_method = "adaptBH",
      gpd_ctrl = gpd_ctrl_custom,
      gamma_ctrl = gamma_ctrl_custom,
      adjust_ctrl = adjust_ctrl_custom
    ), NA
  )
  expect_error(
    compute_p_values(
      obs_stats = obs_mat,
      perm_stats = perm_mat,
      method = "gamma",
      fit_thresh = 0.1,
      alternative = "less",
      null_center = 0,
      adjust_method = "adaptBH",
      gpd_ctrl = gpd_ctrl_custom,
      gamma_ctrl = gamma_ctrl_custom,
      adjust_ctrl = adjust_ctrl_custom
    ), NA
  )
})
