# tests/testthat/test-trajectory.R
#
# Unit tests for gg_trajectory() and its internal helpers
# .parse_initial_conditions() and .build_trajectory_data().


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
std_xlim  <- c(0, 5)
std_ylim  <- c(0, 5)

# Simple linear ODE for integration correctness tests:
# dx/dt = -x, dy/dt = -y -> x(t) = x0*exp(-t), y(t) = y0*exp(-t)
ode_linear_decay <- function(t, y, parameters) {
  list(c(-y[[1L]], -y[[2L]]))
}

# 1D: dy/dt = -y -> y(t) = y0*exp(-t)
ode_decay_1d <- function(t, y, parameters) {
  list(c(-y[[1L]]))
}


# ===========================================================================
# .parse_initial_conditions()
# ===========================================================================

test_that(".parse_initial_conditions() accepts a numeric vector (2D)", {
  result <- ggphasr:::.parse_initial_conditions(c(1, 2), "two.dim")
  expect_type(result, "list")
  expect_length(result, 1L)
  expect_equal(result[[1L]], c(1, 2))
})

test_that(".parse_initial_conditions() accepts a numeric vector (1D)", {
  result <- ggphasr:::.parse_initial_conditions(c(3), "one.dim")
  expect_type(result, "list")
  expect_length(result, 1L)
  expect_equal(result[[1L]], 3)
})

test_that(".parse_initial_conditions() accepts a matrix (2D)", {
  mat    <- matrix(c(1,2, 3,4, 0.5, 1.5), ncol=2, byrow=TRUE)
  result <- ggphasr:::.parse_initial_conditions(mat, "two.dim")
  expect_length(result, 3L)
  expect_equal(result[[1L]], c(1, 2))
  expect_equal(result[[2L]], c(3, 4))
  expect_equal(result[[3L]], c(0.5, 1.5))
})

test_that(".parse_initial_conditions() accepts a list of vectors (2D)", {
  lst    <- list(c(1, 2), c(3, 4))
  result <- ggphasr:::.parse_initial_conditions(lst, "two.dim")
  expect_length(result, 2L)
  expect_equal(result[[1L]], c(1, 2))
  expect_equal(result[[2L]], c(3, 4))
})

test_that(".parse_initial_conditions() errors on wrong vector length (2D)", {
  expect_error(
    ggphasr:::.parse_initial_conditions(c(1, 2, 3), "two.dim"),
    regexp = "length 2"
  )
})

test_that(".parse_initial_conditions() errors on wrong matrix columns", {
  mat <- matrix(1:6, ncol = 3)
  expect_error(
    ggphasr:::.parse_initial_conditions(mat, "two.dim"),
    regexp = "2 column"
  )
})

test_that(".parse_initial_conditions() errors on wrong list element length", {
  lst <- list(c(1, 2), c(1, 2, 3))
  expect_error(
    ggphasr:::.parse_initial_conditions(lst, "two.dim"),
    regexp = "length 2"
  )
})


# ===========================================================================
# .integrate_trajectory()
# ===========================================================================

test_that(".integrate_trajectory() returns a data frame with correct columns (2D)", {
  norm   <- ggphasr:::.normalize_ode(ode_linear_decay, "two.dim")
  result <- ggphasr:::.integrate_trajectory(norm, c(1, 1), 0, 1, 50L,
                                              NULL, "lsoda")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "x", "y") %in% names(result)))
})

test_that(".integrate_trajectory() returns correct columns for 1D", {
  norm   <- ggphasr:::.normalize_ode(ode_decay_1d, "one.dim")
  result <- ggphasr:::.integrate_trajectory(norm, c(2), 0, 1, 50L,
                                              NULL, "lsoda")
  expect_true(all(c("time", "x", "y") %in% names(result)))
})

test_that(".integrate_trajectory() forward: solution matches analytic (2D)", {
  # dx/dt = -x -> x(t) = x0*exp(-t); at t=1, x0=2: x = 2*exp(-1) ≈ 0.7358
  norm   <- ggphasr:::.normalize_ode(ode_linear_decay, "two.dim")
  result <- ggphasr:::.integrate_trajectory(norm, c(2, 3), 0, 1, 200L,
                                              NULL, "lsoda")
  final  <- result[nrow(result), ]
  expect_equal(final$x, 2 * exp(-1), tolerance = 1e-4)
  expect_equal(final$y, 3 * exp(-1), tolerance = 1e-4)
})

test_that(".integrate_trajectory() backward: solution matches analytic", {
  # x(-1) = x0*exp(1) for forward direction;
  # integrating from 0 to -1 gives x(-1) = x0*exp(1)
  norm   <- ggphasr:::.normalize_ode(ode_linear_decay, "two.dim")
  result <- ggphasr:::.integrate_trajectory(norm, c(1, 1), 0, -1, 200L,
                                              NULL, "lsoda")
  # At t = -1: x should be exp(1) ≈ 2.718
  final  <- result[nrow(result), ]
  expect_equal(final$x, exp(1), tolerance = 1e-4)
})

test_that(".integrate_trajectory() returns NULL on integration failure", {
  # A function that always errors
  bad_ode <- function(t, y, parameters) stop("bad ode")
  result  <- ggphasr:::.integrate_trajectory(bad_ode, c(1, 1), 0, 1,
                                               50L, NULL, "lsoda")
  expect_null(result)
})


# ===========================================================================
# gg_trajectory() — return type and structure
# ===========================================================================

test_that("gg_trajectory() returns a list", {
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = c(1, 1),
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  expect_type(result, "list")
})

test_that("gg_trajectory() list elements are ggplot2 layers or scales", {
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = c(1, 1),
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  is_gg <- vapply(result, function(x) {
    inherits(x, "Layer") || inherits(x, "Scale") || inherits(x, "ggproto")
  }, logical(1L))
  expect_true(all(is_gg))
})

test_that("gg_trajectory() composes with gg_flow_field() without error", {
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim       = std_xlim,
                     ylim       = std_ylim,
                     parameters = lv_params) +
    gg_trajectory(ode_lotka_volterra,
                  y0         = c(1, 1),
                  xlim       = std_xlim,
                  ylim       = std_ylim,
                  parameters = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_trajectory() composes with gg_flow_field() + gg_nullclines()", {
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim = std_xlim, ylim = std_ylim,
                     parameters = lv_params) +
    gg_nullclines(ode_lotka_volterra,
                  xlim = std_xlim, ylim = std_ylim,
                  parameters = lv_params) +
    gg_trajectory(ode_lotka_volterra,
                  y0 = c(1, 1), xlim = std_xlim, ylim = std_ylim,
                  parameters = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_trajectory() — multiple initial conditions
# ===========================================================================

test_that("gg_trajectory() accepts a matrix of initial conditions", {
  ics <- matrix(c(0.5,0.5, 1,2, 3,1), ncol=2, byrow=TRUE)
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = ics,
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  expect_type(result, "list")
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim=std_xlim, ylim=std_ylim,
                     parameters=lv_params) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_trajectory() accepts a list of initial conditions", {
  ics <- list(c(0.5, 0.5), c(2, 2), c(3, 1))
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = ics,
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  expect_type(result, "list")
})

test_that("gg_trajectory() uses color scale for multiple ICs", {
  ics    <- list(c(1, 1), c(2, 2))
  result <- gg_trajectory(ode_lotka_volterra,
                           y0 = ics, xlim = std_xlim, ylim = std_ylim,
                           parameters = lv_params)
  # Should include a Scale object for color
  has_scale <- any(vapply(result, function(x) inherits(x, "Scale"),
                          logical(1L)))
  expect_true(has_scale)
})

test_that("gg_trajectory() uses fixed color when color is supplied", {
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = list(c(1,1), c(2,2)),
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params,
                           color      = "steelblue")
  # Should NOT include a color Scale
  has_scale <- any(vapply(result, function(x) inherits(x, "Scale"),
                          logical(1L)))
  expect_false(has_scale)
})


# ===========================================================================
# gg_trajectory() — backward integration
# ===========================================================================

test_that("gg_trajectory() backward integration produces a list", {
  result <- gg_trajectory(ode_example_08,
                           y0           = c(0.1, 2),
                           xlim         = c(-3, 3),
                           ylim         = c(-3, 3),
                           t_end        = 2,
                           t_start_back = -2)
  expect_type(result, "list")
  p <- gg_flow_field(ode_example_08,
                     xlim = c(-3,3), ylim = c(-3,3)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_trajectory() errors when t_start_back is positive", {
  expect_error(
    gg_trajectory(ode_lotka_volterra,
                  y0 = c(1,1), xlim = std_xlim, ylim = std_ylim,
                  parameters = lv_params,
                  t_start_back = 5),
    regexp = "negative"
  )
})


# ===========================================================================
# gg_trajectory() — 1D system
# ===========================================================================

test_that("gg_trajectory() works for a 1D system", {
  result <- gg_trajectory(ode_logistic,
                           y0         = c(1),
                           xlim       = c(0, 10),
                           ylim       = c(-1, 12),
                           system     = "one.dim",
                           parameters = c(r = 1, K = 10))
  expect_type(result, "list")
  p <- gg_flow_field(ode_logistic,
                     xlim = c(0,10), ylim = c(-1,12),
                     system = "one.dim",
                     parameters = c(r=1, K=10)) + result
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_trajectory() — argument options
# ===========================================================================

test_that("gg_trajectory() works with add_arrows = FALSE", {
  result <- gg_trajectory(ode_lotka_volterra,
                           y0         = c(1, 1),
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params,
                           add_arrows = FALSE)
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim=std_xlim, ylim=std_ylim,
                     parameters=lv_params) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_trajectory() works with add_start_point = FALSE", {
  result <- gg_trajectory(ode_lotka_volterra,
                           y0              = c(1, 1),
                           xlim            = std_xlim,
                           ylim            = std_ylim,
                           parameters      = lv_params,
                           add_start_point = FALSE)
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim=std_xlim, ylim=std_ylim,
                     parameters=lv_params) + result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_trajectory() works with Convention B ODE", {
  ode_b2 <- function(x, y, parameters = NULL) c(-x, -y)
  result  <- gg_trajectory(ode_b2,
                            y0   = c(1, 1),
                            xlim = c(-2, 2),
                            ylim = c(-2, 2),
                            t_end = 2)
  expect_type(result, "list")
})


# ===========================================================================
# gg_trajectory() — input validation
# ===========================================================================

test_that("gg_trajectory() errors on invalid xlim", {
  expect_error(
    gg_trajectory(ode_lotka_volterra, y0=c(1,1),
                  xlim=c(5,0), ylim=std_ylim,
                  parameters=lv_params),
    regexp = "xlim"
  )
})

test_that("gg_trajectory() errors on t_end <= 0", {
  expect_error(
    gg_trajectory(ode_lotka_volterra, y0=c(1,1),
                  xlim=std_xlim, ylim=std_ylim,
                  parameters=lv_params, t_end = -1),
    regexp = "t_end"
  )
})

test_that("gg_trajectory() errors on wrong IC dimension", {
  expect_error(
    gg_trajectory(ode_lotka_volterra, y0=c(1,2,3),
                  xlim=std_xlim, ylim=std_ylim,
                  parameters=lv_params),
    regexp = "length 2"
  )
})
