# tests/testthat/test-phase_plane.R
#
# Unit tests for gg_phase_plane() and its internal helpers
# .default_ics(), .equilibrium_shapes(), .equilibrium_fills().


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lv_params   <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
std_xlim    <- c(0, 5)
std_ylim    <- c(0, 5)
log_params  <- c(r = 1, K = 10)


# ===========================================================================
# .default_ics()
# ===========================================================================

test_that(".default_ics() returns a list for 2D", {
  result <- ggphasr:::.default_ics(c(0,5), c(0,5), 3L, "two.dim")
  expect_type(result, "list")
  expect_length(result, 9L)   # 3x3 grid
})

test_that(".default_ics() returns a list for 1D", {
  result <- ggphasr:::.default_ics(c(0,5), c(0,5), 4L, "one.dim")
  expect_type(result, "list")
  expect_length(result, 4L)
})

test_that(".default_ics() each 2D IC has length 2", {
  result <- ggphasr:::.default_ics(c(0,5), c(0,5), 3L, "two.dim")
  lens   <- vapply(result, length, integer(1L))
  expect_true(all(lens == 2L))
})

test_that(".default_ics() each 1D IC has length 1", {
  result <- ggphasr:::.default_ics(c(0,5), c(0,5), 4L, "one.dim")
  lens   <- vapply(result, length, integer(1L))
  expect_true(all(lens == 1L))
})

test_that(".default_ics() ICs are within plot bounds", {
  xlim <- c(0, 5); ylim <- c(-2, 8)
  result <- ggphasr:::.default_ics(xlim, ylim, 4L, "two.dim")
  xs <- vapply(result, function(ic) ic[[1L]], numeric(1L))
  ys <- vapply(result, function(ic) ic[[2L]], numeric(1L))
  expect_true(all(xs >= xlim[[1L]] & xs <= xlim[[2L]]))
  expect_true(all(ys >= ylim[[1L]] & ys <= ylim[[2L]]))
})


# ===========================================================================
# gg_phase_plane() — return type and structure
# ===========================================================================

test_that("gg_phase_plane() returns a ggphasr_result object", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_s3_class(result, "ggphasr_result")
})

test_that("gg_phase_plane() result has $plot and $equilibria components", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_true("plot"       %in% names(result))
  expect_true("equilibria" %in% names(result))
})

test_that("gg_phase_plane() $plot is a ggplot object", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_s3_class(result$plot, "ggplot")
})

test_that("gg_phase_plane() $plot renders without error", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() $equilibria is a data frame when find_equilibria=TRUE", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_s3_class(result$equilibria, "data.frame")
})

test_that("gg_phase_plane() $equilibria is NULL when find_equilibria=FALSE", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            find_equilibria = FALSE)
  expect_null(result$equilibria)
})


# ===========================================================================
# gg_phase_plane() — equilibrium finding
# ===========================================================================

test_that("gg_phase_plane() finds Lotka-Volterra interior equilibrium", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  # Interior equilibrium at (2, 2) should be found
  eq  <- result$equilibria
  has_interior <- any(abs(eq$x - 2) < 0.1 & abs(eq$y - 2) < 0.1)
  expect_true(has_interior)
})

test_that("gg_phase_plane() classifies Lotka-Volterra interior as Center", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  eq <- result$equilibria
  interior_row <- eq[abs(eq$x - 2) < 0.1 & abs(eq$y - 2) < 0.1, ]
  expect_equal(interior_row$classification[[1L]], "Center")
})

test_that("gg_phase_plane() equilibria table has expected columns", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_true(all(c("x", "y", "classification") %in%
                    names(result$equilibria)))
})


# ===========================================================================
# gg_phase_plane() — 1D system
# ===========================================================================

test_that("gg_phase_plane() works for 1D system", {
  result <- gg_phase_plane(ode_logistic,
                            xlim       = c(0, 8),
                            ylim       = c(-1, 12),
                            system     = "one.dim",
                            parameters = log_params)
  expect_s3_class(result, "ggphasr_result")
  expect_s3_class(result$plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() finds logistic equilibria in 1D", {
  result <- gg_phase_plane(ode_logistic,
                            xlim       = c(0, 8),
                            ylim       = c(-1, 12),
                            system     = "one.dim",
                            parameters = log_params)
  eq <- result$equilibria
  expect_true(any(abs(eq$y - 10) < 0.5))   # stable at K=10
  expect_true(any(abs(eq$y - 0)  < 0.5))   # unstable at 0
})


# ===========================================================================
# gg_phase_plane() — component toggles
# ===========================================================================

test_that("gg_phase_plane() works with show_nullclines = FALSE", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            show_nullclines  = FALSE,
                            find_equilibria = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works with show_trajectories = FALSE", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim              = std_xlim,
                            ylim              = std_ylim,
                            parameters        = lv_params,
                            show_trajectories = FALSE,
                            find_equilibria   = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works with custom y0", {
  ics    <- matrix(c(1,1, 2,3, 3,1), ncol=2, byrow=TRUE)
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            y0              = ics,
                            find_equilibria = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works with show_eq_legend = FALSE", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim          = std_xlim,
                            ylim          = std_ylim,
                            parameters    = lv_params,
                            show_eq_legend = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works with backward integration", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            t_start_back    = -5,
                            find_equilibria = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works with proportional arrow_type", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            arrow_type      = "proportional",
                            find_equilibria = FALSE)
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() applies title", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            title           = "Test title",
                            find_equilibria = FALSE)
  built <- ggplot2::ggplot_build(result$plot)
  expect_equal(built$plot$labels$title, "Test title")
})

test_that("gg_phase_plane() applies custom axis labels", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim            = std_xlim,
                            ylim            = std_ylim,
                            parameters      = lv_params,
                            xlab            = "Prey",
                            ylab            = "Predator",
                            find_equilibria = FALSE)
  built <- ggplot2::ggplot_build(result$plot)
  expect_equal(built$plot$labels$x, "Prey")
  expect_equal(built$plot$labels$y, "Predator")
})


# ===========================================================================
# S3 methods
# ===========================================================================

test_that("print.ggphasr_result() runs without error", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  expect_no_error(print(result))
})

test_that("print.ggphasr_result() returns result invisibly", {
  result <- gg_phase_plane(ode_lotka_volterra,
                            xlim       = std_xlim,
                            ylim       = std_ylim,
                            parameters = lv_params)
  ret <- withVisible(print(result))
  expect_false(ret$visible)
})

test_that("add_layer() adds a ggplot2 layer to $plot", {
  result  <- gg_phase_plane(ode_lotka_volterra,
                             xlim            = std_xlim,
                             ylim            = std_ylim,
                             parameters      = lv_params,
                             find_equilibria = FALSE)
  result2 <- add_layer(result, ggplot2::labs(title = "Added title"))
  expect_s3_class(result2, "ggphasr_result")
  built <- ggplot2::ggplot_build(result2$plot)
  expect_equal(built$plot$labels$title, "Added title")
})

test_that("add_layer() preserves $equilibria", {
  result  <- gg_phase_plane(ode_lotka_volterra,
                             xlim       = std_xlim,
                             ylim       = std_ylim,
                             parameters = lv_params)
  result2 <- add_layer(result, ggplot2::labs(title = "Added"))
  expect_s3_class(result2$equilibria, "data.frame")
})

test_that("add_layer() errors on non-ggphasr_result input", {
  expect_error(
    add_layer(list(plot = NULL), ggplot2::labs(title = "x")),
    regexp = "ggphasr_result"
  )
})


# ===========================================================================
# gg_phase_plane() — input validation
# ===========================================================================

test_that("gg_phase_plane() errors on invalid xlim", {
  expect_error(
    gg_phase_plane(ode_lotka_volterra, xlim=c(5,0), ylim=std_ylim,
                   parameters=lv_params),
    regexp = "xlim"
  )
})

test_that("gg_phase_plane() errors on invalid ylim", {
  expect_error(
    gg_phase_plane(ode_lotka_volterra, xlim=std_xlim, ylim=c(5,0),
                   parameters=lv_params),
    regexp = "ylim"
  )
})


# ===========================================================================
# gg_phase_plane() — additional systems
# ===========================================================================

test_that("gg_phase_plane() works for Van der Pol", {
  result <- gg_phase_plane(ode_van_der_pol,
                            xlim            = c(-3, 3),
                            ylim            = c(-4, 4),
                            parameters      = c(mu = 1),
                            t_end           = 15,
                            find_equilibria = TRUE)
  expect_s3_class(result, "ggphasr_result")
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works for competition model", {
  result <- gg_phase_plane(ode_competition,
                            xlim       = c(0, 12),
                            ylim       = c(0, 12),
                            parameters = c(r1=1,r2=1,K1=10,K2=10,
                                           a12=0.5,a21=0.5))
  expect_s3_class(result, "ggphasr_result")
  expect_no_error(ggplot2::ggplot_build(result$plot))
})

test_that("gg_phase_plane() works for Convention B ODE", {
  ode_b <- function(x, y, parameters = NULL) c(-x, -y)
  result <- gg_phase_plane(ode_b,
                            xlim            = c(-2, 2),
                            ylim            = c(-2, 2),
                            find_equilibria = FALSE)
  expect_s3_class(result, "ggphasr_result")
  expect_no_error(ggplot2::ggplot_build(result$plot))
})
