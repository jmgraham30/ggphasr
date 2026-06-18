# tests/testthat/test-nullclines.R
#
# Unit tests for gg_nullclines() and its internal helper
# .compute_nullcline_grid().


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lv_params  <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
std_xlim   <- c(0, 5)
std_ylim   <- c(0, 5)

# A system with analytically known nullclines:
# ode_example_07: dx/dt = -x, dy/dt = -y
# x-nullcline: x = 0 (the y-axis)
# y-nullcline: y = 0 (the x-axis)
ode_easy_nullclines <- function(t, y, parameters) {
  list(c(-y[[1L]], -y[[2L]]))
}


# ===========================================================================
# .compute_nullcline_grid()
# ===========================================================================

test_that(".compute_nullcline_grid() returns correct columns for 2D", {
  norm <- ggphasr:::.normalize_ode(ode_easy_nullclines, "two.dim")
  df   <- ggphasr:::.compute_nullcline_grid(norm, "two.dim",
                                             c(-2, 2), c(-2, 2),
                                             20L, NULL)
  expect_true(all(c("x", "y", "f", "g") %in% names(df)))
  expect_equal(nrow(df), 400L)
})

test_that(".compute_nullcline_grid() returns correct columns for 1D", {
  norm <- ggphasr:::.normalize_ode(ode_example_01, "one.dim")
  df   <- ggphasr:::.compute_nullcline_grid(norm, "one.dim",
                                             c(0, 4), c(-4, 4),
                                             50L, NULL)
  expect_true(all(c("y", "f") %in% names(df)))
  expect_equal(nrow(df), 50L)
})

test_that(".compute_nullcline_grid() f column is correct for simple system", {
  # ode_easy_nullclines: dx/dt = -x, dy/dt = -y
  # f should equal -x at each grid point
  norm <- ggphasr:::.normalize_ode(ode_easy_nullclines, "two.dim")
  df   <- ggphasr:::.compute_nullcline_grid(norm, "two.dim",
                                             c(-2, 2), c(-2, 2),
                                             5L, NULL)
  expect_equal(df$f, -df$x, tolerance = 1e-12)
  expect_equal(df$g, -df$y, tolerance = 1e-12)
})

test_that(".compute_nullcline_grid() f = 0 where nullcline crosses for 1D", {
  # ode_example_01: dy/dt = 4 - y^2 = 0 at y = +-2
  norm <- ggphasr:::.normalize_ode(ode_example_01, "one.dim")
  df   <- ggphasr:::.compute_nullcline_grid(norm, "one.dim",
                                             c(0, 4), c(-3, 3),
                                             200L, NULL)
  # Find rows closest to y = 2 and y = -2
  idx_pos <- which.min(abs(df$y - 2))
  idx_neg <- which.min(abs(df$y + 2))
  expect_equal(df$f[[idx_pos]], 0, tolerance = 0.1)
  expect_equal(df$f[[idx_neg]], 0, tolerance = 0.1)
})


# ===========================================================================
# gg_nullclines() — return type and structure
# ===========================================================================

test_that("gg_nullclines() returns a list", {
  result <- gg_nullclines(ode_lotka_volterra,
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  expect_type(result, "list")
})

test_that("gg_nullclines() list elements are ggplot2 layers or scales", {
  result <- gg_nullclines(ode_lotka_volterra,
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params)
  is_gg <- vapply(result, function(x) {
    inherits(x, "Layer") || inherits(x, "Scale") || inherits(x, "ggproto")
  }, logical(1L))
  expect_true(all(is_gg))
})

test_that("gg_nullclines() composes with gg_flow_field() without error", {
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim       = std_xlim,
                     ylim       = std_ylim,
                     parameters = lv_params) +
    gg_nullclines(ode_lotka_volterra,
                  xlim       = std_xlim,
                  ylim       = std_ylim,
                  parameters = lv_params)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_nullclines() works with add_legend = FALSE", {
  result <- gg_nullclines(ode_lotka_volterra,
                           xlim       = std_xlim,
                           ylim       = std_ylim,
                           parameters = lv_params,
                           add_legend = FALSE)
  expect_type(result, "list")
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim = std_xlim, ylim = std_ylim,
                     parameters = lv_params) +
    result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_nullclines() works for 1D system", {
  result <- gg_nullclines(ode_logistic,
                           xlim       = c(0, 4),
                           ylim       = c(-2, 12),
                           system     = "one.dim",
                           parameters = c(r = 1, K = 10))
  expect_type(result, "list")
  p <- gg_flow_field(ode_logistic,
                     xlim = c(0,4), ylim = c(-2, 12),
                     system = "one.dim",
                     parameters = c(r = 1, K = 10)) +
    result
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_nullclines() accepts Convention B ODE", {
  ode_b2 <- function(x, y, parameters = NULL) c(y, -x)
  result  <- gg_nullclines(ode_b2,
                            xlim = c(-2, 2), ylim = c(-2, 2))
  expect_type(result, "list")
})


# ===========================================================================
# gg_nullclines() — input validation
# ===========================================================================

test_that("gg_nullclines() errors on invalid xlim", {
  expect_error(
    gg_nullclines(ode_lotka_volterra, xlim = c(5, 0), ylim = std_ylim,
                  parameters = lv_params),
    regexp = "xlim"
  )
})

test_that("gg_nullclines() errors on invalid ylim", {
  expect_error(
    gg_nullclines(ode_lotka_volterra, xlim = std_xlim, ylim = c(5, 0),
                  parameters = lv_params),
    regexp = "ylim"
  )
})


# ===========================================================================
# gg_nullclines() — visual correctness
# ===========================================================================

test_that("gg_nullclines() produces layers that render for competition model", {
  # Competition model has linear nullclines — a good rendering test
  comp_params <- c(r1=1, r2=1, K1=10, K2=10, a12=0.5, a21=0.5)
  p <- gg_flow_field(ode_competition,
                     xlim = c(0,15), ylim = c(0,15),
                     parameters = comp_params) +
    gg_nullclines(ode_competition,
                  xlim = c(0,15), ylim = c(0,15),
                  parameters = comp_params)
  expect_no_error(ggplot2::ggplot_build(p))
})
