# tests/testthat/test-ode_systems_1d.R
#
# Unit tests for the four built-in 1D ODE systems.
#
# Each function is tested for:
#   (1) Correct return structure (list containing a length-1 numeric vector)
#   (2) Correct derivative values at analytically known points
#   (3) Correct behavior at equilibrium points (derivative = 0)
#   (4) Named and positional parameter passing both work
#   (5) Default parameters produce a valid result


# ---------------------------------------------------------------------------
# Helper: check that output has the right structure
# ---------------------------------------------------------------------------

expect_ode_structure <- function(result, dim = 1L) {
  expect_type(result, "list")
  expect_length(result, 1L)
  expect_type(result[[1L]], "double")
  expect_length(result[[1L]], dim)
}


# ===========================================================================
# ode_exponential()
# ===========================================================================

test_that("ode_exponential() returns correct structure", {
  result <- ode_exponential(t = 0, y = c(1), parameters = c(r = 0.5))
  expect_ode_structure(result)
})

test_that("ode_exponential() computes dy/dt = r * y correctly", {
  # r = 2, y = 3 -> dy/dt = 6
  result <- ode_exponential(t = 0, y = c(3), parameters = c(r = 2))
  expect_equal(result[[1L]], 6)
})

test_that("ode_exponential() equilibrium at y = 0", {
  result <- ode_exponential(t = 0, y = c(0), parameters = c(r = 1))
  expect_equal(result[[1L]], 0)
})

test_that("ode_exponential() models decay when r < 0", {
  # r = -1, y = 5 -> dy/dt = -5
  result <- ode_exponential(t = 0, y = c(5), parameters = c(r = -1))
  expect_equal(result[[1L]], -5)
})

test_that("ode_exponential() works with positional parameters", {
  # parameters[[1]] treated as r
  result <- ode_exponential(t = 0, y = c(2), parameters = c(0.5))
  expect_equal(result[[1L]], 1)
})

test_that("ode_exponential() default parameters produce valid output", {
  result <- ode_exponential(t = 0, y = c(1))
  expect_ode_structure(result)
  # default r = 0.5, y = 1 -> dy/dt = 0.5
  expect_equal(result[[1L]], 0.5)
})

test_that("ode_exponential() is time-invariant (autonomous)", {
  r1 <- ode_exponential(t = 0,   y = c(2), parameters = c(r = 1))
  r2 <- ode_exponential(t = 100, y = c(2), parameters = c(r = 1))
  expect_equal(r1[[1L]], r2[[1L]])
})


# ===========================================================================
# ode_logistic()
# ===========================================================================

test_that("ode_logistic() returns correct structure", {
  result <- ode_logistic(t = 0, y = c(5), parameters = c(r = 1, K = 10))
  expect_ode_structure(result)
})

test_that("ode_logistic() computes dy/dt = r*y*(1 - y/K) correctly", {
  # r = 1, K = 10, y = 5 -> dy/dt = 1 * 5 * (1 - 0.5) = 2.5
  result <- ode_logistic(t = 0, y = c(5), parameters = c(r = 1, K = 10))
  expect_equal(result[[1L]], 2.5)
})

test_that("ode_logistic() unstable equilibrium at y = 0", {
  result <- ode_logistic(t = 0, y = c(0), parameters = c(r = 1, K = 10))
  expect_equal(result[[1L]], 0)
})

test_that("ode_logistic() stable equilibrium at y = K", {
  result <- ode_logistic(t = 0, y = c(10), parameters = c(r = 1, K = 10))
  expect_equal(result[[1L]], 0)
})

test_that("ode_logistic() gives negative derivative above K (overshoot)", {
  # y > K -> population declines back toward K
  result <- ode_logistic(t = 0, y = c(12), parameters = c(r = 1, K = 10))
  expect_lt(result[[1L]], 0)
})

test_that("ode_logistic() gives maximum growth rate near y = K/2", {
  # Growth rate is maximized at y = K/2 for standard logistic
  # dy/dt at K/2 should exceed dy/dt at K/4 and dy/dt at 3K/4
  r1 <- ode_logistic(t = 0, y = c(5),  parameters = c(r = 1, K = 10))[[1L]]
  r2 <- ode_logistic(t = 0, y = c(2.5),parameters = c(r = 1, K = 10))[[1L]]
  r3 <- ode_logistic(t = 0, y = c(7.5),parameters = c(r = 1, K = 10))[[1L]]
  expect_gt(r1, r2)
  expect_gt(r1, r3)
})

test_that("ode_logistic() works with positional parameters", {
  result_named    <- ode_logistic(t = 0, y = c(5), parameters = c(r = 2, K = 8))
  result_positional <- ode_logistic(t = 0, y = c(5), parameters = c(2, 8))
  expect_equal(result_named[[1L]], result_positional[[1L]])
})

test_that("ode_logistic() default parameters produce valid output", {
  result <- ode_logistic(t = 0, y = c(5))
  expect_ode_structure(result)
  # default r = 1, K = 10, y = 5 -> dy/dt = 2.5
  expect_equal(result[[1L]], 2.5)
})


# ===========================================================================
# ode_monomolecular()
# ===========================================================================

test_that("ode_monomolecular() returns correct structure", {
  result <- ode_monomolecular(t = 0, y = c(4), parameters = c(r = 1, K = 10))
  expect_ode_structure(result)
})

test_that("ode_monomolecular() computes dy/dt = r*(K - y) correctly", {
  # r = 1, K = 10, y = 4 -> dy/dt = 1 * (10 - 4) = 6
  result <- ode_monomolecular(t = 0, y = c(4), parameters = c(r = 1, K = 10))
  expect_equal(result[[1L]], 6)
})

test_that("ode_monomolecular() stable equilibrium at y = K", {
  result <- ode_monomolecular(t = 0, y = c(10), parameters = c(r = 1, K = 10))
  expect_equal(result[[1L]], 0)
})

test_that("ode_monomolecular() has monotonically decreasing growth rate", {
  # Growth rate decreases as y increases toward K
  r1 <- ode_monomolecular(t = 0, y = c(1),  parameters = c(r = 1, K = 10))[[1L]]
  r2 <- ode_monomolecular(t = 0, y = c(5),  parameters = c(r = 1, K = 10))[[1L]]
  r3 <- ode_monomolecular(t = 0, y = c(9),  parameters = c(r = 1, K = 10))[[1L]]
  expect_gt(r1, r2)
  expect_gt(r2, r3)
})

test_that("ode_monomolecular() gives negative derivative above K", {
  result <- ode_monomolecular(t = 0, y = c(12), parameters = c(r = 1, K = 10))
  expect_lt(result[[1L]], 0)
})

test_that("ode_monomolecular() works with positional parameters", {
  result_named      <- ode_monomolecular(t = 0, y = c(3), parameters = c(r = 2, K = 8))
  result_positional <- ode_monomolecular(t = 0, y = c(3), parameters = c(2, 8))
  expect_equal(result_named[[1L]], result_positional[[1L]])
})

test_that("ode_monomolecular() default parameters produce valid output", {
  result <- ode_monomolecular(t = 0, y = c(4))
  expect_ode_structure(result)
  # default r = 1, K = 10, y = 4 -> dy/dt = 6
  expect_equal(result[[1L]], 6)
})


# ===========================================================================
# ode_von_bertalanffy()
# ===========================================================================

test_that("ode_von_bertalanffy() returns correct structure", {
  result <- ode_von_bertalanffy(t = 0, y = c(8), parameters = c(alpha = 1, beta = 0.5))
  expect_ode_structure(result)
})

test_that("ode_von_bertalanffy() computes dy/dt = alpha*y^(2/3) - beta*y correctly", {
  # alpha = 1, beta = 0.5, y = 1 -> dy/dt = 1*1^(2/3) - 0.5*1 = 0.5
  result <- ode_von_bertalanffy(t = 0, y = c(1), parameters = c(alpha = 1, beta = 0.5))
  expect_equal(result[[1L]], 0.5)
})

test_that("ode_von_bertalanffy() equilibrium at y = (alpha/beta)^3", {
  # alpha = 1, beta = 0.5 -> y* = (1/0.5)^3 = 8
  # dy/dt = 1*8^(2/3) - 0.5*8 = 4 - 4 = 0
  result <- ode_von_bertalanffy(t = 0, y = c(8), parameters = c(alpha = 1, beta = 0.5))
  expect_equal(result[[1L]], 0, tolerance = 1e-10)
})

test_that("ode_von_bertalanffy() equilibrium at y = 0", {
  result <- ode_von_bertalanffy(t = 0, y = c(0), parameters = c(alpha = 1, beta = 0.5))
  expect_equal(result[[1L]], 0)
})

test_that("ode_von_bertalanffy() gives positive growth below equilibrium", {
  # y < y* = 8 -> dy/dt > 0
  result <- ode_von_bertalanffy(t = 0, y = c(4), parameters = c(alpha = 1, beta = 0.5))
  expect_gt(result[[1L]], 0)
})

test_that("ode_von_bertalanffy() gives negative growth above equilibrium", {
  # y > y* = 8 -> dy/dt < 0
  result <- ode_von_bertalanffy(t = 0, y = c(12), parameters = c(alpha = 1, beta = 0.5))
  expect_lt(result[[1L]], 0)
})

test_that("ode_von_bertalanffy() works with positional parameters", {
  result_named      <- ode_von_bertalanffy(t = 0, y = c(1),
                                            parameters = c(alpha = 2, beta = 1))
  result_positional <- ode_von_bertalanffy(t = 0, y = c(1),
                                            parameters = c(2, 1))
  expect_equal(result_named[[1L]], result_positional[[1L]])
})

test_that("ode_von_bertalanffy() default parameters produce valid output", {
  result <- ode_von_bertalanffy(t = 0, y = c(1))
  expect_ode_structure(result)
  # default alpha = 1, beta = 0.5, y = 1 -> dy/dt = 0.5
  expect_equal(result[[1L]], 0.5)
})

test_that("ode_von_bertalanffy() equilibrium location scales correctly with parameters", {
  # y* = (alpha/beta)^3, so doubling alpha while keeping beta fixed
  # should move the equilibrium to (2/0.5)^3 = 64
  result <- ode_von_bertalanffy(t = 0, y = c(64),
                                 parameters = c(alpha = 2, beta = 0.5))
  expect_equal(result[[1L]], 0, tolerance = 1e-10)
})


# ===========================================================================
# Cross-function tests
# ===========================================================================

test_that("all 1D ODE functions are compatible with .normalize_ode()", {
  fns <- list(
    ode_exponential     = list(fn = ode_exponential,     params = c(r = 1)),
    ode_logistic        = list(fn = ode_logistic,        params = c(r = 1, K = 10)),
    ode_monomolecular   = list(fn = ode_monomolecular,   params = c(r = 1, K = 10)),
    ode_von_bertalanffy = list(fn = ode_von_bertalanffy, params = c(alpha = 1, beta = 0.5))
  )

  for (name in names(fns)) {
    fn     <- fns[[name]]$fn
    params <- fns[[name]]$params
    normalized <- ggphasr:::.normalize_ode(fn, system = "one.dim")
    result <- normalized(t = 0, y = c(1), parameters = params)
    expect_type(result, "list")
    expect_length(result[[1L]], 1L)
  }
})

test_that("all 1D ODE functions pass .validate_ode()", {
  fns <- list(
    list(fn = ode_exponential,     params = c(r = 1)),
    list(fn = ode_logistic,        params = c(r = 1, K = 10)),
    list(fn = ode_monomolecular,   params = c(r = 1, K = 10)),
    list(fn = ode_von_bertalanffy, params = c(alpha = 1, beta = 0.5))
  )

  for (item in fns) {
    normalized <- ggphasr:::.normalize_ode(item$fn, system = "one.dim")
    expect_true(
      ggphasr:::.validate_ode(normalized,
                               system     = "one.dim",
                               parameters = item$params)
    )
  }
})
