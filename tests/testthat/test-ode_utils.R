# tests/testthat/test-ode_utils.R
#
# Unit tests for the internal ODE utility functions in ode_utils.R.
#
# These tests use simple, analytically understood ODE systems so that
# expected outputs can be computed by hand and do not depend on any other
# ggphasr functions.

# ---------------------------------------------------------------------------
# Helper ODE functions used across multiple tests
# ---------------------------------------------------------------------------

# Convention A — 1D: dy/dt = r*y  (exponential growth)
ode_a_1d <- function(t, y, parameters) {
  r <- parameters[1L]
  list(c(r * y[1L]))
}

# Convention A — 2D: Lotka-Volterra
ode_a_2d <- function(t, y, parameters) {
  x <- y[1L]
  v <- y[2L]
  a <- parameters[1L]; b <- parameters[2L]
  c <- parameters[3L]; d <- parameters[4L]
  list(c(
    a * x - b * x * v,
    c * x * v - d * v
  ))
}

# Convention B1 — 1D: dy/dt = -y  (exponential decay)
ode_b1 <- function(y, parameters = NULL) {
  -y
}

# Convention B2 — 2D: dx/dt = y, dy/dt = -x  (harmonic oscillator)
ode_b2 <- function(x, y, parameters = NULL) {
  c(y, -x)
}

# A non-function (to test error handling)
not_a_function <- 42L

# A function with no arguments (to test error handling)
ode_no_args <- function() 0

# A function that returns a numeric vector instead of a list (bad Convention A)
ode_bad_return <- function(t, y, parameters) {
  c(y[1L], y[2L])   # missing list() wrapper
}

# A function that errors during evaluation
ode_errors <- function(t, y, parameters) {
  stop("intentional error in ODE")
}

# A Convention A 2D function that returns wrong-length output
ode_wrong_length <- function(t, y, parameters) {
  list(c(1, 2, 3))  # 3 derivatives for a 2D system
}

# ---------------------------------------------------------------------------
# .detect_ode_convention()
# ---------------------------------------------------------------------------

test_that(".detect_ode_convention() correctly identifies Convention A", {
  expect_equal(ggphasr:::.detect_ode_convention(ode_a_1d), "A")
  expect_equal(ggphasr:::.detect_ode_convention(ode_a_2d), "A")
})

test_that(".detect_ode_convention() correctly identifies Convention B1 (1D)", {
  expect_equal(ggphasr:::.detect_ode_convention(ode_b1), "B1")
})

test_that(".detect_ode_convention() correctly identifies Convention B2 (2D)", {
  expect_equal(ggphasr:::.detect_ode_convention(ode_b2), "B2")
})

test_that(".detect_ode_convention() errors if `deriv` is not a function", {
  expect_error(
    ggphasr:::.detect_ode_convention(not_a_function),
    regexp = "must be a function"
  )
})

test_that(".detect_ode_convention() errors if `deriv` has no arguments", {
  expect_error(
    ggphasr:::.detect_ode_convention(ode_no_args),
    regexp = "at least one argument"
  )
})

test_that(".detect_ode_convention() errors on unrecognized argument names", {
  ode_weird <- function(state, params) list(c(0, 0))
  expect_error(
    ggphasr:::.detect_ode_convention(ode_weird),
    regexp = "Could not detect"
  )
})

# ---------------------------------------------------------------------------
# .normalize_ode()
# ---------------------------------------------------------------------------

test_that(".normalize_ode() returns Convention A functions unchanged", {
  normalized <- ggphasr:::.normalize_ode(ode_a_2d, system = "two.dim")
  # Should be identical to the original function object
  expect_identical(normalized, ode_a_2d)
})

test_that(".normalize_ode() wraps Convention B1 into a valid Convention A function", {
  normalized <- ggphasr:::.normalize_ode(ode_b1, system = "one.dim")
  result <- normalized(t = 0, y = c(2.0), parameters = NULL)

  expect_type(result, "list")
  expect_length(result[[1L]], 1L)
  expect_equal(result[[1L]], -2.0)  # ode_b1: dy/dt = -y, y = 2 -> -2
})

test_that(".normalize_ode() wraps Convention B2 into a valid Convention A function", {
  normalized <- ggphasr:::.normalize_ode(ode_b2, system = "two.dim")
  result <- normalized(t = 0, y = c(1.0, 0.0), parameters = NULL)

  expect_type(result, "list")
  expect_length(result[[1L]], 2L)
  # ode_b2: dx/dt = y = 0, dy/dt = -x = -1
  expect_equal(result[[1L]], c(0.0, -1.0))
})

test_that(".normalize_ode() errors when B2 function used with one.dim", {
  expect_error(
    ggphasr:::.normalize_ode(ode_b2, system = "one.dim"),
    regexp = "2D Convention B"
  )
})

test_that(".normalize_ode() errors when B1 function used with two.dim", {
  expect_error(
    ggphasr:::.normalize_ode(ode_b1, system = "two.dim"),
    regexp = "1D Convention B"
  )
})

# ---------------------------------------------------------------------------
# .validate_ode()
# ---------------------------------------------------------------------------

test_that(".validate_ode() passes silently for a valid 2D Convention A ODE", {
  normalized <- ggphasr:::.normalize_ode(ode_a_2d, system = "two.dim")
  expect_true(
    ggphasr:::.validate_ode(normalized, system = "two.dim",
                             parameters = c(1, 0.5, 0.5, 1))
  )
})

test_that(".validate_ode() passes silently for a valid 1D Convention A ODE", {
  normalized <- ggphasr:::.normalize_ode(ode_a_1d, system = "one.dim")
  expect_true(
    ggphasr:::.validate_ode(normalized, system = "one.dim",
                             parameters = c(0.5))
  )
})

test_that(".validate_ode() errors when `deriv` throws during evaluation", {
  normalized <- ggphasr:::.normalize_ode(ode_errors, system = "two.dim")
  expect_error(
    ggphasr:::.validate_ode(normalized, system = "two.dim"),
    regexp = "threw an error"
  )
})

test_that(".validate_ode() errors when `deriv` returns non-list", {
  normalized <- ggphasr:::.normalize_ode(ode_bad_return, system = "two.dim")
  expect_error(
    ggphasr:::.validate_ode(normalized, system = "two.dim"),
    regexp = "must return a list"
  )
})

test_that(".validate_ode() errors when derivative vector has wrong length", {
  normalized <- ggphasr:::.normalize_ode(ode_wrong_length, system = "two.dim")
  expect_error(
    ggphasr:::.validate_ode(normalized, system = "two.dim"),
    regexp = "length 2"
  )
})

# ---------------------------------------------------------------------------
# .eval_ode()
# ---------------------------------------------------------------------------

test_that(".eval_ode() returns a plain numeric vector, not a list", {
  normalized <- ggphasr:::.normalize_ode(ode_a_2d, system = "two.dim")
  result <- ggphasr:::.eval_ode(normalized,
                                 t          = 0,
                                 y          = c(2.0, 1.0),
                                 parameters = c(1, 0.5, 0.5, 1))
  expect_type(result, "double")
  expect_false(is.list(result))
})

test_that(".eval_ode() returns correct values for Lotka-Volterra", {
  # At (x=2, y=1) with (a=1, b=0.5, c=0.5, d=1):
  # dx/dt = 1*2 - 0.5*2*1 = 2 - 1 = 1
  # dy/dt = 0.5*2*1 - 1*1 = 1 - 1 = 0
  normalized <- ggphasr:::.normalize_ode(ode_a_2d, system = "two.dim")
  result <- ggphasr:::.eval_ode(normalized,
                                 t          = 0,
                                 y          = c(2.0, 1.0),
                                 parameters = c(1, 0.5, 0.5, 1))
  expect_equal(result, c(1.0, 0.0))
})

test_that(".eval_ode() returns correct values for exponential growth", {
  # dy/dt = r*y, at y=3, r=2 -> dy/dt = 6
  normalized <- ggphasr:::.normalize_ode(ode_a_1d, system = "one.dim")
  result <- ggphasr:::.eval_ode(normalized,
                                 t          = 0,
                                 y          = c(3.0),
                                 parameters = c(2.0))
  expect_equal(result, 6.0)
})

test_that(".eval_ode() works correctly for wrapped Convention B2 functions", {
  # Harmonic oscillator: dx/dt = y, dy/dt = -x
  # At (x=0, y=1): dx/dt = 1, dy/dt = 0
  normalized <- ggphasr:::.normalize_ode(ode_b2, system = "two.dim")
  result <- ggphasr:::.eval_ode(normalized,
                                 t          = 0,
                                 y          = c(0.0, 1.0),
                                 parameters = NULL)
  expect_equal(result, c(1.0, 0.0))
})


# ===========================================================================
# .legend_theme()
# ===========================================================================

test_that(".legend_theme() returns a theme object for standard positions", {
  for (pos in c("right", "left", "top", "bottom", "none")) {
    th <- ggphasr:::.legend_theme(pos)
    expect_s3_class(th, "theme")
    expect_equal(th$legend.position, pos)
  }
})

test_that('.legend_theme("inside") returns a theme with inside position', {
  th <- ggphasr:::.legend_theme("inside")
  expect_s3_class(th, "theme")
  expect_equal(th$legend.position, "inside")
  expect_true(!is.null(th$legend.position.inside))
})

test_that(".legend_theme() accepts numeric c(x,y) vector", {
  th <- ggphasr:::.legend_theme(c(0.9, 0.9))
  expect_s3_class(th, "theme")
  expect_equal(th$legend.position, "inside")
  expect_equal(th$legend.position.inside, c(0.9, 0.9))
})

test_that(".legend_theme() can be added to a ggplot", {
  p <- ggplot2::ggplot() + ggphasr:::.legend_theme("none")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
