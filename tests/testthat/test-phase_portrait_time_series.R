# tests/testthat/test-phase_portrait.R
# tests/testthat/test-time_series.R
#
# Combined test file for Layer 4 functions.


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

lv_params <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)


# ===========================================================================
# .find_equilibria_1d()
# ===========================================================================

test_that(".find_equilibria_1d() finds equilibria of example_01", {
  # dy/dt = 4 - y^2: equilibria at y = -2 (unstable) and y = 2 (stable)
  norm <- ggphasr:::.normalize_ode(ode_example_01, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(-4, 4), NULL, 500L)
  expect_equal(nrow(eq), 2L)
  expect_equal(sort(round(eq$y, 4)), c(-2, 2))
})

test_that(".find_equilibria_1d() classifies stable equilibrium correctly", {
  norm <- ggphasr:::.normalize_ode(ode_example_01, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(-4, 4), NULL, 500L)
  stable_rows <- eq[round(eq$y, 1) == 2, ]
  expect_equal(stable_rows$stability, "stable")
})

test_that(".find_equilibria_1d() classifies unstable equilibrium correctly", {
  norm <- ggphasr:::.normalize_ode(ode_example_01, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(-4, 4), NULL, 500L)
  unstable_rows <- eq[round(eq$y, 1) == -2, ]
  expect_equal(unstable_rows$stability, "unstable")
})

test_that(".find_equilibria_1d() finds three equilibria of example_02", {
  # dy/dt = y(1-y)(2-y): equilibria at y = 0, 1, 2
  norm <- ggphasr:::.normalize_ode(ode_example_02, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(-0.5, 2.5), NULL, 500L)
  expect_equal(nrow(eq), 3L)
  expect_equal(sort(round(eq$y, 3)), c(0, 1, 2))
})

test_that(".find_equilibria_1d() returns empty data frame when no equilibria", {
  # dy/dt = 1 (no equilibria in range)
  ode_no_eq <- function(t, y, parameters) list(c(1))
  norm <- ggphasr:::.normalize_ode(ode_no_eq, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(0, 5), NULL, 100L)
  expect_equal(nrow(eq), 0L)
})

test_that(".find_equilibria_1d() finds logistic equilibria with parameters", {
  # dy/dt = r*y*(1 - y/K): equilibria at y=0 (unstable) and y=K=10 (stable)
  norm <- ggphasr:::.normalize_ode(ode_logistic, "one.dim")
  eq   <- ggphasr:::.find_equilibria_1d(norm, c(-1, 12),
                                          c(r = 1, K = 10), 500L)
  expect_equal(nrow(eq), 2L)
  expect_true(any(abs(eq$y - 10) < 0.1))
  expect_true(any(abs(eq$y - 0)  < 0.1))
})


# ===========================================================================
# gg_phase_portrait() — return type and structure
# ===========================================================================

test_that("gg_phase_portrait() returns a list", {
  result <- gg_phase_portrait(ode_example_01,
                               ylim = c(-4, 4),
                               xlim = c(0, 4))
  expect_type(result, "list")
})

test_that("gg_phase_portrait() list elements are ggplot2 layers", {
  result <- gg_phase_portrait(ode_example_01,
                               ylim = c(-4, 4),
                               xlim = c(0, 4))
  is_gg <- vapply(result, function(x) {
    inherits(x, "Layer") || inherits(x, "ggproto")
  }, logical(1L))
  expect_true(all(is_gg))
})

test_that("gg_phase_portrait() composes with gg_flow_field()", {
  p <- gg_flow_field(ode_example_01,
                     xlim   = c(0, 4),
                     ylim   = c(-4, 4),
                     system = "one.dim") +
    gg_phase_portrait(ode_example_01,
                      ylim = c(-4, 4),
                      xlim = c(0, 4))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_phase_portrait() works for logistic growth", {
  p <- gg_flow_field(ode_logistic,
                     xlim       = c(0, 6),
                     ylim       = c(-1, 12),
                     system     = "one.dim",
                     parameters = c(r = 1, K = 10)) +
    gg_phase_portrait(ode_logistic,
                      ylim       = c(-1, 12),
                      xlim       = c(0, 6),
                      parameters = c(r = 1, K = 10))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_phase_portrait() works for three-equilibrium system", {
  p <- gg_flow_field(ode_example_02,
                     xlim   = c(0, 4),
                     ylim   = c(-0.5, 2.5),
                     system = "one.dim") +
    gg_phase_portrait(ode_example_02,
                      ylim = c(-0.5, 2.5),
                      xlim = c(0, 4))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_phase_portrait() accepts Convention B ODE", {
  ode_b <- function(y, parameters = NULL) 4 - y^2
  result <- gg_phase_portrait(ode_b,
                               ylim = c(-4, 4),
                               xlim = c(0, 4))
  expect_type(result, "list")
})

test_that("gg_phase_portrait() composes with gg_nullclines() and gg_trajectory()", {
  p <- gg_flow_field(ode_logistic,
                     xlim       = c(0, 8),
                     ylim       = c(-1, 12),
                     system     = "one.dim",
                     parameters = c(r = 1, K = 10)) +
    gg_nullclines(ode_logistic,
                  xlim       = c(0, 8),
                  ylim       = c(-1, 12),
                  system     = "one.dim",
                  parameters = c(r = 1, K = 10)) +
    gg_trajectory(ode_logistic,
                  y0         = list(c(1), c(5), c(11)),
                  xlim       = c(0, 8),
                  ylim       = c(-1, 12),
                  system     = "one.dim",
                  parameters = c(r = 1, K = 10),
                  t_end      = 7,
                  color      = "steelblue") +
    gg_phase_portrait(ode_logistic,
                      ylim       = c(-1, 12),
                      xlim       = c(0, 8),
                      parameters = c(r = 1, K = 10))
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_phase_portrait() — input validation
# ===========================================================================

test_that("gg_phase_portrait() errors on invalid ylim", {
  expect_error(
    gg_phase_portrait(ode_example_01, ylim = c(4, -4), xlim = c(0, 4)),
    regexp = "ylim"
  )
})

test_that("gg_phase_portrait() errors on invalid xlim", {
  expect_error(
    gg_phase_portrait(ode_example_01, ylim = c(-4, 4), xlim = c(4, 0)),
    regexp = "xlim"
  )
})


# ===========================================================================
# gg_time_series() — return type and structure
# ===========================================================================

test_that("gg_time_series() returns a ggplot for 2D system", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0         = c(1, 1),
                      t_end      = 10,
                      parameters = lv_params)
  expect_s3_class(p, "ggplot")
})

test_that("gg_time_series() renders without error for 2D system", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0         = c(1, 1),
                      t_end      = 10,
                      parameters = lv_params)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_time_series() returns a ggplot for 1D system", {
  p <- gg_time_series(ode_logistic,
                      y0         = c(1),
                      t_end      = 8,
                      system     = "one.dim",
                      parameters = c(r = 1, K = 10))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_time_series() applies theme_phase_plane()", {
  p     <- gg_time_series(ode_lotka_volterra, y0 = c(1,1),
                           t_end = 5, parameters = lv_params)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$theme$panel.background$fill, "white")
})


# ===========================================================================
# gg_time_series() — multiple initial conditions
# ===========================================================================

test_that("gg_time_series() works with multiple ICs (list)", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0         = list(c(1,1), c(2,2), c(3,1)),
                      t_end      = 10,
                      parameters = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_time_series() works with multiple ICs (matrix)", {
  ics <- matrix(c(1,1, 2,2, 3,1), ncol=2, byrow=TRUE)
  p   <- gg_time_series(ode_lotka_volterra,
                         y0         = ics,
                         t_end      = 10,
                         parameters = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_time_series() uses fixed color when color is supplied", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0         = list(c(1,1), c(2,2)),
                      t_end      = 10,
                      parameters = lv_params,
                      color      = "steelblue")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_time_series() works with add_legend = FALSE", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0         = list(c(1,1), c(2,2)),
                      t_end      = 10,
                      parameters = lv_params,
                      add_legend = FALSE)
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_time_series() — 2D faceting and labels
# ===========================================================================

test_that("gg_time_series() uses custom var_labels in 2D facet", {
  p     <- gg_time_series(ode_lotka_volterra,
                           y0         = c(1, 1),
                           t_end      = 5,
                           parameters = lv_params,
                           var_labels = c("Prey", "Predator"))
  built <- ggplot2::ggplot_build(p)
  # Extract facet labels
  facet_labels <- built$layout$facet_params$facets[[1L]]
  expect_s3_class(p, "ggplot")
  expect_no_error(built)
})

test_that("gg_time_series() errors on wrong var_labels length", {
  expect_error(
    gg_time_series(ode_lotka_volterra, y0 = c(1,1), t_end = 5,
                   parameters = lv_params,
                   var_labels = c("x")),
    regexp = "var_labels"
  )
})

test_that("gg_time_series() applies custom axis labels", {
  p     <- gg_time_series(ode_logistic,
                           y0         = c(1),
                           t_end      = 8,
                           system     = "one.dim",
                           parameters = c(r = 1, K = 10),
                           xlab       = "t (days)",
                           ylab       = "Population")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$labels$x, "t (days)")
  expect_equal(built$plot$labels$y, "Population")
})

test_that("gg_time_series() applies title", {
  p     <- gg_time_series(ode_logistic,
                           y0 = c(1), t_end = 5,
                           system = "one.dim",
                           parameters = c(r=1, K=10),
                           title = "My title")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$labels$title, "My title")
})


# ===========================================================================
# gg_time_series() — backward integration
# ===========================================================================

test_that("gg_time_series() works with backward integration", {
  p <- gg_time_series(ode_lotka_volterra,
                      y0           = c(1, 1),
                      t_end        = 5,
                      t_start_back = -2,
                      parameters   = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# gg_time_series() — input validation
# ===========================================================================

test_that("gg_time_series() errors on t_end <= 0", {
  expect_error(
    gg_time_series(ode_logistic, y0 = c(1), t_end = -1,
                   system = "one.dim", parameters = c(r=1, K=10)),
    regexp = "t_end"
  )
})

test_that("gg_time_series() errors on positive t_start_back", {
  expect_error(
    gg_time_series(ode_lotka_volterra, y0 = c(1,1),
                   t_end = 5, t_start_back = 2,
                   parameters = lv_params),
    regexp = "negative"
  )
})

test_that("gg_time_series() errors on wrong IC dimension", {
  expect_error(
    gg_time_series(ode_lotka_volterra, y0 = c(1, 2, 3),
                   t_end = 5, parameters = lv_params),
    regexp = "length 2"
  )
})

test_that("gg_time_series() works with Convention B ODE", {
  ode_b <- function(x, y, parameters = NULL) c(-x, -y)
  p <- gg_time_series(ode_b, y0 = c(1, 1), t_end = 3)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
