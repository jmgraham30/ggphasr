# tests/testthat/test-manifolds.R
#
# Unit tests for gg_manifolds() and its internal helpers
# .get_manifold_eigenvectors() and .integrate_manifold_branch().


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)

# example_08: dx/dt = x, dy/dt = -y
# Saddle at origin; unstable eigenvector = (1,0), stable eigenvector = (0,1)
# Analytic Jacobian: [[1,0],[0,-1]]

# example_11: nonlinear competition, saddle at (1,1)
eq11    <- c(1, 1)
eq11_cl <- classify_equilibrium(ode_example_11, equilibrium = eq11)


# ===========================================================================
# .get_manifold_eigenvectors()
# ===========================================================================

test_that(".get_manifold_eigenvectors() returns expected list elements", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   NULL, 1e-6)
  expect_true(all(c("stable_vec", "unstable_vec",
                    "lambda_s", "lambda_u", "J") %in% names(result)))
})

test_that(".get_manifold_eigenvectors() eigenvectors are unit vectors", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   NULL, 1e-6)
  expect_equal(sqrt(sum(result$stable_vec^2)),   1, tolerance = 1e-10)
  expect_equal(sqrt(sum(result$unstable_vec^2)), 1, tolerance = 1e-10)
})

test_that(".get_manifold_eigenvectors() stable eigenvalue is negative", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   NULL, 1e-6)
  expect_lt(Re(result$lambda_s), 0)
})

test_that(".get_manifold_eigenvectors() unstable eigenvalue is positive", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   NULL, 1e-6)
  expect_gt(Re(result$lambda_u), 0)
})

test_that(".get_manifold_eigenvectors() correct eigenvalues for example_08", {
  # Jacobian [[1,0],[0,-1]]: eigenvalues 1 and -1
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   NULL, 1e-6)
  expect_equal(Re(result$lambda_u),  1, tolerance = 1e-5)
  expect_equal(Re(result$lambda_s), -1, tolerance = 1e-5)
})

test_that(".get_manifold_eigenvectors() uses pre-computed Jacobian correctly", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  eq_cl  <- classify_equilibrium(ode_example_08, equilibrium = c(0,0))
  result <- ggphasr:::.get_manifold_eigenvectors(norm, c(0,0), NULL,
                                                   eq_cl, 1e-6)
  # Should still give eigenvalues 1 and -1
  expect_equal(Re(result$lambda_u),  1, tolerance = 1e-5)
  expect_equal(Re(result$lambda_s), -1, tolerance = 1e-5)
})


# ===========================================================================
# .integrate_manifold_branch()
# ===========================================================================

test_that(".integrate_manifold_branch() returns a data frame with x, y", {
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  seed   <- c(0.1, 0)
  result <- ggphasr:::.integrate_manifold_branch(norm, seed, 2, 100L,
                                                   NULL, "lsoda")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("x", "y") %in% names(result)))
})

test_that(".integrate_manifold_branch() forward integration moves away correctly", {
  # example_08: along unstable direction (1,0), x grows exponentially
  # seed at (0.1, 0), forward t=1: x should be ~0.1*exp(1) = 0.272
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.integrate_manifold_branch(norm, c(0.1, 0),
                                                   1, 200L, NULL, "lsoda")
  expect_equal(tail(result$x, 1), 0.1 * exp(1), tolerance = 1e-3)
  expect_equal(tail(result$y, 1), 0,             tolerance = 1e-3)
})

test_that(".integrate_manifold_branch() backward integration converges to saddle", {
  # example_08: along stable direction (0,1), backward integration
  # from (0, 0.1) at t=-1: y should approach 0.1*exp(1) = 0.272
  # (going backward on -y means y increases toward 0 from the saddle perspective)
  norm   <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  result <- ggphasr:::.integrate_manifold_branch(norm, c(0, 0.1),
                                                   -1, 200L, NULL, "lsoda")
  expect_equal(tail(result$y, 1), 0.1 * exp(1), tolerance = 1e-3)
})

test_that(".integrate_manifold_branch() returns NULL on failure", {
  bad_ode <- function(t, y, parameters) stop("error")
  result  <- ggphasr:::.integrate_manifold_branch(bad_ode, c(0.1, 0),
                                                    2, 50L, NULL, "lsoda")
  expect_null(result)
})


# ===========================================================================
# gg_manifolds() — return type and structure
# ===========================================================================

test_that("gg_manifolds() returns a list", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0, 0))
  expect_type(result, "list")
})

test_that("gg_manifolds() list elements are ggplot2 layers or scales", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0, 0))
  is_gg  <- vapply(result, function(x) {
    inherits(x, "Layer") || inherits(x, "Scale") || inherits(x, "ggproto") || inherits(x, "theme")
  }, logical(1L))
  expect_true(all(is_gg))
})

test_that("gg_manifolds() composes with gg_flow_field() without error", {
  p <- gg_flow_field(ode_example_08,
                     xlim = c(-3, 3), ylim = c(-3, 3)) +
    gg_manifolds(ode_example_08, equilibrium = c(0, 0))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() renders without error for example_11 saddle", {
  p <- gg_flow_field(ode_example_11,
                     xlim = c(0, 4), ylim = c(0, 3)) +
    gg_manifolds(ode_example_11,
                 equilibrium   = eq11,
                 eq_classified = eq11_cl,
                 t_manifold    = 4)
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_manifolds() — argument options
# ===========================================================================

test_that("gg_manifolds() works with draw_stable = FALSE", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0,0),
                          draw_stable = FALSE)
  expect_type(result, "list")
  p <- gg_flow_field(ode_example_08, xlim=c(-3,3), ylim=c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() works with draw_unstable = FALSE", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0,0),
                          draw_unstable = FALSE)
  expect_type(result, "list")
  p <- gg_flow_field(ode_example_08, xlim=c(-3,3), ylim=c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() works with add_legend = FALSE", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0,0),
                          add_legend = FALSE)
  p <- gg_flow_field(ode_example_08, xlim=c(-3,3), ylim=c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() works with add_arrows = FALSE", {
  result <- gg_manifolds(ode_example_08, equilibrium = c(0,0),
                          add_arrows = FALSE)
  p <- gg_flow_field(ode_example_08, xlim=c(-3,3), ylim=c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() works with pre-computed eq_classified", {
  eq_cl  <- classify_equilibrium(ode_example_08, equilibrium = c(0,0))
  result <- gg_manifolds(ode_example_08, equilibrium = c(0,0),
                          eq_classified = eq_cl)
  expect_type(result, "list")
  p <- gg_flow_field(ode_example_08, xlim=c(-3,3), ylim=c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_manifolds() accepts Convention B ODE", {
  ode_b <- function(x, y, parameters = NULL) c(x, -y)
  result <- gg_manifolds(ode_b, equilibrium = c(0, 0))
  expect_type(result, "list")
})


# ===========================================================================
# gg_manifolds() — non-saddle warning
# ===========================================================================

test_that("gg_manifolds() warns and returns empty list for non-saddle (with eq_classified)", {
  # example_07 origin is a stable node, not a saddle
  eq_cl  <- classify_equilibrium(ode_example_07, equilibrium = c(0,0))
  expect_warning(
    result <- gg_manifolds(ode_example_07, equilibrium = c(0,0),
                            eq_classified = eq_cl),
    regexp = "not.*Saddle|Saddle"
  )
  expect_length(result, 0L)
})

test_that("gg_manifolds() warns and returns empty list for non-saddle (no eq_classified)", {
  # example_09 origin is a stable spiral
  expect_warning(
    result <- gg_manifolds(ode_example_09, equilibrium = c(0,0)),
    regexp = "saddle|Saddle"
  )
  expect_length(result, 0L)
})


# ===========================================================================
# gg_manifolds() — input validation
# ===========================================================================

test_that("gg_manifolds() errors on equilibrium of wrong length", {
  expect_error(
    gg_manifolds(ode_example_08, equilibrium = c(0, 0, 0)),
    regexp = "length 2"
  )
})

test_that("gg_manifolds() errors when both draw_stable and draw_unstable are FALSE", {
  expect_error(
    gg_manifolds(ode_example_08, equilibrium = c(0,0),
                 draw_stable = FALSE, draw_unstable = FALSE),
    regexp = "At least one"
  )
})

test_that("gg_manifolds() errors on non-positive t_manifold", {
  expect_error(
    gg_manifolds(ode_example_08, equilibrium = c(0,0), t_manifold = -1),
    regexp = "t_manifold"
  )
})

test_that("gg_manifolds() errors on invalid eq_classified", {
  expect_error(
    gg_manifolds(ode_example_08, equilibrium = c(0,0),
                 eq_classified = data.frame(x = 1)),
    regexp = "classify_equilibrium"
  )
})


# ===========================================================================
# Full workflow test
# ===========================================================================

test_that("find + classify + manifolds full workflow renders correctly", {
  eq_list <- find_equilibrium(ode_example_11,
                               y0 = c(0.8, 0.8))
  eq_cl   <- classify_equilibrium(ode_example_11,
                                   equilibrium = eq_list[[1L]])

  p <- gg_flow_field(ode_example_11,
                     xlim = c(0, 4), ylim = c(0, 3)) +
    gg_nullclines(ode_example_11,
                  xlim = c(0, 4), ylim = c(0, 3)) +
    gg_manifolds(ode_example_11,
                 equilibrium   = eq_list[[1L]],
                 eq_classified = eq_cl,
                 t_manifold    = 4)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
