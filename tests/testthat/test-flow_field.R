# tests/testthat/test-flow_field.R
#
# Unit tests for gg_flow_field().
#
# Strategy:
#   - Return type and class checks (every path through the function)
#   - Input validation (bad xlim, ylim, n_points)
#   - Data correctness: verify .compute_flow_field() and .scale_arrows()
#     produce the right structure and values
#   - All major argument combinations: system, arrow_type, color options,
#     multi-system, origin lines, title
#   - Composability: gg_flow_field() output accepts + operator


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Simple Convention A 2D ODE for tests (linear saddle, easy to reason about)
ode_linear_saddle <- function(t, y, parameters) {
  list(c(y[[1L]], -y[[2L]]))
}

# Simple Convention A 1D ODE for tests
ode_linear_1d <- function(t, y, parameters) {
  list(c(-y[[1L]]))
}

# Convention B 2D ODE (harmonic oscillator) for convention detection tests
ode_harmonic_b <- function(x, y, parameters = NULL) {
  c(y, -x)
}

std_params_lv <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)


# ===========================================================================
# Return type
# ===========================================================================

test_that("gg_flow_field() returns a ggplot object for 2D system", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim = c(-2, 2), ylim = c(-2, 2))
  expect_s3_class(p, "ggplot")
})

test_that("gg_flow_field() returns a ggplot object for 1D system", {
  p <- gg_flow_field(ode_linear_1d,
                     xlim   = c(0, 4),
                     ylim   = c(-3, 3),
                     system = "one.dim")
  expect_s3_class(p, "ggplot")
})

test_that("gg_flow_field() output renders without error", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim = c(-2, 2), ylim = c(-2, 2))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() accepts Convention B ODE functions", {
  p <- gg_flow_field(ode_harmonic_b,
                     xlim = c(-2, 2), ylim = c(-2, 2))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# Input validation
# ===========================================================================

test_that("gg_flow_field() errors on invalid xlim", {
  expect_error(
    gg_flow_field(ode_linear_saddle, xlim = c(2, -2), ylim = c(-2, 2)),
    regexp = "xlim"
  )
  expect_error(
    gg_flow_field(ode_linear_saddle, xlim = c(1, 1), ylim = c(-2, 2)),
    regexp = "xlim"
  )
})

test_that("gg_flow_field() errors on invalid ylim", {
  expect_error(
    gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(2, -2)),
    regexp = "ylim"
  )
})

test_that("gg_flow_field() errors on n_points < 2", {
  expect_error(
    gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(-2, 2),
                  n_points = 1),
    regexp = "n_points"
  )
})

test_that("gg_flow_field() errors when list deriv has unnamed elements", {
  expect_error(
    gg_flow_field(list(ode_linear_saddle, ode_harmonic_b),
                  xlim = c(-2, 2), ylim = c(-2, 2)),
    regexp = "named"
  )
})


# ===========================================================================
# .compute_flow_field() internals
# ===========================================================================

test_that(".compute_flow_field() returns correct columns for 2D", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       xlim = c(-1, 1), ylim = c(-1, 1),
                                       n_points = 5L,
                                       parameters = NULL)
  expect_true(all(c("x", "y", "dx", "dy") %in% names(df)))
  expect_equal(nrow(df), 25L)   # 5 * 5 grid
})

test_that(".compute_flow_field() returns correct columns for 1D", {
  norm <- ggphasr:::.normalize_ode(ode_linear_1d, system = "one.dim")
  df <- ggphasr:::.compute_flow_field(norm, "one.dim",
                                       xlim = c(0, 4), ylim = c(-2, 2),
                                       n_points = 5L,
                                       parameters = NULL)
  expect_true(all(c("x", "y", "dx", "dy") %in% names(df)))
  expect_equal(nrow(df), 25L)
  # For 1D, dx should always be 1
  expect_true(all(df$dx == 1))
})

test_that(".compute_flow_field() gives correct derivative at known 2D point", {
  # ode_linear_saddle: dx/dt = x, dy/dt = -y
  # At (1, 1): dx = 1, dy = -1
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       xlim = c(0, 2), ylim = c(0, 2),
                                       n_points = 3L,
                                       parameters = NULL)
  row <- df[df$x == 1 & df$y == 1, ]
  expect_equal(row$dx, 1)
  expect_equal(row$dy, -1)
})

test_that(".compute_flow_field() gives correct derivative for 1D system", {
  # ode_linear_1d: dy/dt = -y; at y=2, dy = -2
  norm <- ggphasr:::.normalize_ode(ode_linear_1d, system = "one.dim")
  df <- ggphasr:::.compute_flow_field(norm, "one.dim",
                                       xlim = c(0, 4), ylim = c(0, 4),
                                       n_points = 3L,
                                       parameters = NULL)
  row <- df[abs(df$y - 2) < 1e-10, ][1L, ]
  expect_equal(row$dy, -2)
})


# ===========================================================================
# .scale_arrows() internals
# ===========================================================================

test_that(".scale_arrows() adds xend and yend columns", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       c(-1,1), c(-1,1), 5L, NULL)
  df2 <- ggphasr:::.scale_arrows(df, "equal", c(-1,1), c(-1,1), 5L)
  expect_true(all(c("xend", "yend", "magnitude") %in% names(df2)))
})

test_that(".scale_arrows() equal: all non-zero arrows have same length", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       c(-1,1), c(-1,1), 5L, NULL)
  df2 <- ggphasr:::.scale_arrows(df, "equal", c(-1,1), c(-1,1), 5L)
  lengths <- sqrt(df2$dx_scaled^2 + df2$dy_scaled^2)
  # Exclude equilibrium at origin (magnitude = 0)
  non_zero <- df2$magnitude > .Machine$double.eps
  expect_equal(length(unique(round(lengths[non_zero], 10))), 1L)
})

test_that(".scale_arrows() proportional: longer arrows have larger magnitude", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       c(0.5, 2), c(0.5, 2), 5L, NULL)
  df2 <- ggphasr:::.scale_arrows(df, "proportional",
                                   c(0.5, 2), c(0.5, 2), 5L)
  display_len <- sqrt(df2$dx_scaled^2 + df2$dy_scaled^2)
  # Rank correlation between magnitude and display length should be 1
  expect_gt(cor(df2$magnitude, display_len, method = "spearman"), 0.99)
})

test_that(".scale_arrows() proportional with max_magnitude scales correctly", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       c(0.5, 2), c(0.5, 2), 3L, NULL)
  df_auto  <- ggphasr:::.scale_arrows(df, "proportional",
                                       c(0.5, 2), c(0.5, 2), 3L)
  df_fixed <- ggphasr:::.scale_arrows(df, "proportional",
                                       c(0.5, 2), c(0.5, 2), 3L,
                                       max_magnitude = 10)
  # With a larger fixed max, arrows should be shorter
  mean_auto  <- mean(sqrt(df_auto$dx_scaled^2  + df_auto$dy_scaled^2))
  mean_fixed <- mean(sqrt(df_fixed$dx_scaled^2 + df_fixed$dy_scaled^2))
  expect_gt(mean_auto, mean_fixed)
})

test_that(".scale_arrows() arrows never exceed cell size", {
  norm <- ggphasr:::.normalize_ode(ode_linear_saddle, system = "two.dim")
  df <- ggphasr:::.compute_flow_field(norm, "two.dim",
                                       c(0.5, 2), c(0.5, 2), 5L, NULL)
  df2 <- ggphasr:::.scale_arrows(df, "equal", c(0.5, 2), c(0.5, 2), 5L,
                                   arrow_length_scale = 0.85)
  cell_size <- min(diff(c(0.5, 2)), diff(c(0.5, 2))) / (5 - 1)
  max_len   <- cell_size * 0.85
  display_len <- sqrt(df2$dx_scaled^2 + df2$dy_scaled^2)
  expect_true(all(display_len <= max_len + 1e-10))
})


# ===========================================================================
# Argument combinations
# ===========================================================================

test_that("gg_flow_field() works with arrow_type = 'proportional'", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim       = c(-2, 2),
                     ylim       = c(-2, 2),
                     arrow_type = "proportional")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() works with color_by_magnitude = TRUE", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim              = c(-2, 2),
                     ylim              = c(-2, 2),
                     color_by_magnitude = TRUE)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() works with add_origin_lines = FALSE", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim            = c(-2, 2),
                     ylim            = c(-2, 2),
                     add_origin_lines = FALSE)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() works with a title", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim  = c(-2, 2),
                     ylim  = c(-2, 2),
                     title = "Test plot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$labels$title, "Test plot")
})

test_that("gg_flow_field() applies default x-axis label correctly", {
  p_2d <- gg_flow_field(ode_linear_saddle,
                         xlim = c(-2, 2), ylim = c(-2, 2))
  p_1d <- gg_flow_field(ode_linear_1d,
                         xlim   = c(0, 4), ylim = c(-2, 2),
                         system = "one.dim")
  expect_equal(ggplot2::ggplot_build(p_2d)$plot$labels$x, "x")
  expect_equal(ggplot2::ggplot_build(p_1d)$plot$labels$x, "t")
})

test_that("gg_flow_field() applies custom axis labels", {
  p <- gg_flow_field(ode_lotka_volterra,
                     xlim       = c(0, 5),
                     ylim       = c(0, 5),
                     parameters = std_params_lv,
                     xlab       = "Prey",
                     ylab       = "Predator")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$labels$x, "Prey")
  expect_equal(built$plot$labels$y, "Predator")
})

test_that("gg_flow_field() works with n_points = 5", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim     = c(-2, 2),
                     ylim     = c(-2, 2),
                     n_points = 5L)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() works with max_magnitude supplied", {
  p <- gg_flow_field(ode_linear_saddle,
                     xlim          = c(-2, 2),
                     ylim          = c(-2, 2),
                     arrow_type    = "proportional",
                     max_magnitude = 5)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# Multi-system
# ===========================================================================

test_that("gg_flow_field() works with multiple systems (named list)", {
  p <- gg_flow_field(
    deriv = list(
      saddle    = ode_linear_saddle,
      harmonic  = ode_example_06
    ),
    xlim = c(-2, 2),
    ylim = c(-2, 2)
  )
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() multi-system adds a legend", {
  p <- gg_flow_field(
    deriv = list(saddle = ode_linear_saddle, harmonic = ode_example_06),
    xlim = c(-2, 2), ylim = c(-2, 2)
  )
  built <- ggplot2::ggplot_build(p)
  # A legend should be present (color aesthetic mapped)
  expect_true("colour" %in% names(built$plot$mapping) ||
              any(vapply(built$plot$layers, function(l) {
                "colour" %in% names(l$mapping)
              }, logical(1L))))
})

test_that("gg_flow_field() multi-system with per-system parameters", {
  p <- gg_flow_field(
    deriv = list(lv1 = ode_lotka_volterra, lv2 = ode_lotka_volterra),
    xlim = c(0, 5), ylim = c(0, 5),
    parameters = list(
      lv1 = c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1),
      lv2 = c(alpha = 2, beta = 0.5, delta = 0.5, gamma = 1)
    )
  )
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# ===========================================================================
# Composability
# ===========================================================================

test_that("gg_flow_field() output accepts + ggplot2::labs()", {
  p <- gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::labs(title = "Added title", x = "State x")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() output accepts + ggplot2::theme() overrides", {
  p <- gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "lightyellow"))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() output accepts + geom_point() for equilibria", {
  p <- gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(-2, 2)) +
    ggplot2::geom_point(data = data.frame(x = 0, y = 0),
                        ggplot2::aes(x = x, y = y),
                        size = 3, color = "red")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("gg_flow_field() applies theme_phase_plane() by default", {
  p <- gg_flow_field(ode_linear_saddle, xlim = c(-2, 2), ylim = c(-2, 2))
  built <- ggplot2::ggplot_build(p)
  # theme_phase_plane sets white panel background
  bg_fill <- built$plot$theme$panel.background$fill
  expect_equal(bg_fill, "white")
})
