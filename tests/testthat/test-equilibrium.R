# tests/testthat/test-equilibrium.R
#
# Unit tests for find_equilibrium() and classify_equilibrium(), plus
# internal helpers .numerical_jacobian(), .classify_2d(), .classify_1d(),
# and .deduplicate_equilibria().
#
# All expected equilibrium locations and classifications are analytically
# verified from the ODE definitions.


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lv_params   <- c(alpha = 1, beta = 0.5, delta = 0.5, gamma = 1)
comp_params <- c(r1 = 1, r2 = 1, K1 = 10, K2 = 10, a12 = 0.5, a21 = 0.5)

# Tolerance for equilibrium location comparisons
eq_tol <- 1e-4


# ===========================================================================
# .numerical_jacobian()
# ===========================================================================

test_that(".numerical_jacobian() returns a matrix of correct dimensions (2D)", {
  norm <- ggphasr:::.normalize_ode(ode_lotka_volterra, "two.dim")
  J    <- ggphasr:::.numerical_jacobian(norm, c(2, 2), lv_params)
  expect_true(is.matrix(J))
  expect_equal(dim(J), c(2L, 2L))
})

test_that(".numerical_jacobian() returns a 1x1 matrix for 1D system", {
  norm <- ggphasr:::.normalize_ode(ode_logistic, "one.dim")
  J    <- ggphasr:::.numerical_jacobian(norm, c(5), c(r = 1, K = 10))
  expect_true(is.matrix(J))
  expect_equal(dim(J), c(1L, 1L))
})

test_that(".numerical_jacobian() is correct for linear system (example_07)", {
  # dx/dt = -x, dy/dt = -y -> J = [[-1, 0], [0, -1]] everywhere
  norm <- ggphasr:::.normalize_ode(ode_example_07, "two.dim")
  J    <- ggphasr:::.numerical_jacobian(norm, c(0, 0), NULL)
  expect_equal(J, matrix(c(-1, 0, 0, -1), nrow = 2L), tolerance = 1e-5)
})

test_that(".numerical_jacobian() is correct for saddle (example_08)", {
  # dx/dt = x, dy/dt = -y -> J = [[1, 0], [0, -1]] everywhere
  norm <- ggphasr:::.normalize_ode(ode_example_08, "two.dim")
  J    <- ggphasr:::.numerical_jacobian(norm, c(0, 0), NULL)
  expect_equal(J, matrix(c(1, 0, 0, -1), nrow = 2L), tolerance = 1e-5)
})

test_that(".numerical_jacobian() is correct for 1D logistic at equilibrium", {
  # dy/dt = r*y*(1 - y/K); at y=K: df/dy = r*(1 - 2*y/K) = r*(1-2) = -r
  norm <- ggphasr:::.normalize_ode(ode_logistic, "one.dim")
  J    <- ggphasr:::.numerical_jacobian(norm, c(10), c(r = 1, K = 10))
  expect_equal(J[[1L, 1L]], -1, tolerance = 1e-5)
})


# ===========================================================================
# .classify_2d()
# ===========================================================================

test_that(".classify_2d() identifies stable node", {
  # Eigenvalues: -1, -2 (both real, negative)
  expect_equal(ggphasr:::.classify_2d(c(-1, -2)), "Stable node")
})

test_that(".classify_2d() identifies unstable node", {
  expect_equal(ggphasr:::.classify_2d(c(1, 2)), "Unstable node")
})

test_that(".classify_2d() identifies saddle", {
  expect_equal(ggphasr:::.classify_2d(c(-1, 2)), "Saddle")
})

test_that(".classify_2d() identifies stable spiral", {
  expect_equal(ggphasr:::.classify_2d(c(-1 + 2i, -1 - 2i)), "Stable spiral")
})

test_that(".classify_2d() identifies unstable spiral", {
  expect_equal(ggphasr:::.classify_2d(c(1 + 2i, 1 - 2i)), "Unstable spiral")
})

test_that(".classify_2d() identifies center", {
  expect_equal(ggphasr:::.classify_2d(c(0 + 2i, 0 - 2i)), "Center")
})

test_that(".classify_2d() identifies non-isolated equilibrium (det = 0)", {
  expect_equal(ggphasr:::.classify_2d(c(0, -1)), "Non-isolated equilibrium")
})


# ===========================================================================
# .classify_1d()
# ===========================================================================

test_that(".classify_1d() identifies stable equilibrium (df/dy < 0)", {
  expect_equal(ggphasr:::.classify_1d(-0.5), "Stable")
})

test_that(".classify_1d() identifies unstable equilibrium (df/dy > 0)", {
  expect_equal(ggphasr:::.classify_1d(0.5), "Unstable")
})

test_that(".classify_1d() identifies inconclusive case (df/dy = 0)", {
  expect_equal(ggphasr:::.classify_1d(0), "Inconclusive (df/dy = 0)")
})


# ===========================================================================
# .deduplicate_equilibria()
# ===========================================================================

test_that(".deduplicate_equilibria() keeps distinct roots", {
  roots  <- list(c(0, 0), c(1, 0), c(0, 1))
  result <- ggphasr:::.deduplicate_equilibria(roots, tol = 1e-4)
  expect_length(result, 3L)
})

test_that(".deduplicate_equilibria() removes near-duplicate roots", {
  roots  <- list(c(0, 0), c(1e-6, 1e-6), c(1, 0))
  result <- ggphasr:::.deduplicate_equilibria(roots, tol = 1e-4)
  expect_length(result, 2L)
})

test_that(".deduplicate_equilibria() handles single-element list", {
  roots  <- list(c(1, 2))
  result <- ggphasr:::.deduplicate_equilibria(roots)
  expect_length(result, 1L)
})


# ===========================================================================
# find_equilibrium() — single initial guess
# ===========================================================================

test_that("find_equilibrium() returns a list", {
  result <- find_equilibrium(ode_lotka_volterra,
                              y0         = c(1.5, 1.5),
                              parameters = lv_params)
  expect_type(result, "list")
})

test_that("find_equilibrium() finds Lotka-Volterra interior equilibrium", {
  # Interior equilibrium at (gamma/delta, alpha/beta) = (2, 2)
  result <- find_equilibrium(ode_lotka_volterra,
                              y0         = c(1.5, 1.5),
                              parameters = lv_params)
  expect_length(result, 1L)
  expect_equal(result[[1L]], c(2, 2), tolerance = eq_tol)
})

test_that("find_equilibrium() finds logistic equilibrium at K (1D)", {
  result <- find_equilibrium(ode_logistic,
                              y0         = c(8),
                              system     = "one.dim",
                              parameters = c(r = 1, K = 10))
  expect_length(result, 1L)
  expect_equal(result[[1L]][[1L]], 10, tolerance = eq_tol)
})

test_that("find_equilibrium() finds trivial equilibrium at origin", {
  result <- find_equilibrium(ode_lotka_volterra,
                              y0         = c(0.01, 0.01),
                              parameters = lv_params)
  expect_length(result, 1L)
  expect_equal(result[[1L]], c(0, 0), tolerance = eq_tol)
})

test_that("find_equilibrium() finds example_07 equilibrium at origin", {
  result <- find_equilibrium(ode_example_07, y0 = c(1, 1))
  expect_length(result, 1L)
  expect_equal(result[[1L]], c(0, 0), tolerance = eq_tol)
})


# ===========================================================================
# find_equilibrium() — grid search (y0 = NULL)
# ===========================================================================

test_that("find_equilibrium() grid search finds all equilibria of example_11", {
  # example_11 has four equilibria: (0,0), (3,0), (0,2), (1,1)
  result <- find_equilibrium(ode_example_11,
                              y0     = NULL,
                              xlim   = c(0, 4),
                              ylim   = c(0, 3),
                              n_grid = 12L)
  expect_gte(length(result), 4L)
  # Check that (1,1) is among the results
  found_coexist <- any(vapply(result, function(r) {
    sqrt(sum((r - c(1, 1))^2)) < eq_tol * 10
  }, logical(1L)))
  expect_true(found_coexist)
})

test_that("find_equilibrium() grid search works for 1D system", {
  # Logistic: equilibria at y=0 and y=K=10
  result <- find_equilibrium(ode_logistic,
                              y0         = NULL,
                              system     = "one.dim",
                              ylim       = c(-1, 12),
                              parameters = c(r = 1, K = 10),
                              n_grid     = 15L)
  expect_gte(length(result), 2L)
  ys <- sort(vapply(result, function(r) r[[1L]], numeric(1L)))
  expect_equal(ys[[1L]], 0,  tolerance = eq_tol)
  expect_equal(ys[[length(ys)]], 10, tolerance = eq_tol)
})

test_that("find_equilibrium() errors when y0=NULL and ylim missing", {
  expect_error(
    find_equilibrium(ode_logistic, y0 = NULL, system = "one.dim"),
    regexp = "ylim"
  )
})

test_that("find_equilibrium() errors when y0=NULL and xlim missing for 2D", {
  expect_error(
    find_equilibrium(ode_lotka_volterra, y0 = NULL, ylim = c(0, 5)),
    regexp = "xlim"
  )
})

test_that("find_equilibrium() accepts matrix of starting points", {
  starts <- matrix(c(0.1, 0.1, 1.5, 1.5), ncol = 2, byrow = TRUE)
  result <- find_equilibrium(ode_lotka_volterra,
                              y0         = starts,
                              parameters = lv_params)
  expect_type(result, "list")
  expect_gte(length(result), 1L)
})


# ===========================================================================
# classify_equilibrium() — return structure
# ===========================================================================

test_that("classify_equilibrium() returns a data frame", {
  result <- classify_equilibrium(ode_lotka_volterra,
                                  equilibrium = c(2, 2),
                                  parameters  = lv_params)
  expect_s3_class(result, "data.frame")
})

test_that("classify_equilibrium() returns exactly one row", {
  result <- classify_equilibrium(ode_lotka_volterra,
                                  equilibrium = c(2, 2),
                                  parameters  = lv_params)
  expect_equal(nrow(result), 1L)
})

test_that("classify_equilibrium() has expected columns for 2D system", {
  result <- classify_equilibrium(ode_lotka_volterra,
                                  equilibrium = c(2, 2),
                                  parameters  = lv_params)
  expected_cols <- c("x", "y", "classification", "tr", "det",
                      "jacobian_11", "jacobian_12", "jacobian_21",
                      "jacobian_22", "lambda_1_re", "lambda_1_im",
                      "lambda_2_re", "lambda_2_im")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("classify_equilibrium() has expected columns for 1D system", {
  result <- classify_equilibrium(ode_logistic,
                                  equilibrium = c(10),
                                  system      = "one.dim",
                                  parameters  = c(r = 1, K = 10))
  expect_true("classification" %in% names(result))
  expect_true("jacobian_11"    %in% names(result))
  expect_true(is.na(result$x))
  expect_true(is.na(result$tr))
})


# ===========================================================================
# classify_equilibrium() — correctness
# ===========================================================================

test_that("classify_equilibrium() identifies Lotka-Volterra center", {
  result <- classify_equilibrium(ode_lotka_volterra,
                                  equilibrium = c(2, 2),
                                  parameters  = lv_params)
  expect_equal(result$classification, "Center")
  expect_equal(result$x, 2, tolerance = 1e-6)
  expect_equal(result$y, 2, tolerance = 1e-6)
})

test_that("classify_equilibrium() identifies stable node (example_07)", {
  # dx/dt = -x, dy/dt = -y: eigenvalues -1, -1 -> stable node
  result <- classify_equilibrium(ode_example_07, equilibrium = c(0, 0))
  expect_equal(result$classification, "Stable node")
})

test_that("classify_equilibrium() identifies saddle (example_08)", {
  # dx/dt = x, dy/dt = -y: eigenvalues 1, -1 -> saddle
  result <- classify_equilibrium(ode_example_08, equilibrium = c(0, 0))
  expect_equal(result$classification, "Saddle")
})

test_that("classify_equilibrium() identifies stable spiral (example_09)", {
  # dx/dt = -x+y, dy/dt = -x-y: eigenvalues -1±i -> stable spiral
  result <- classify_equilibrium(ode_example_09, equilibrium = c(0, 0))
  expect_equal(result$classification, "Stable spiral")
})

test_that("classify_equilibrium() identifies saddle (example_10)", {
  # dx/dt = x+y, dy/dt = x-y: eigenvalues ±sqrt(2) -> saddle
  result <- classify_equilibrium(ode_example_10, equilibrium = c(0, 0))
  expect_equal(result$classification, "Saddle")
})

test_that("classify_equilibrium() identifies stable logistic equilibrium (1D)", {
  # At y=K=10: df/dy = r*(1-2) = -r < 0 -> stable
  result <- classify_equilibrium(ode_logistic,
                                  equilibrium = c(10),
                                  system      = "one.dim",
                                  parameters  = c(r = 1, K = 10))
  expect_equal(result$classification, "Stable")
  expect_equal(result$y, 10)
})

test_that("classify_equilibrium() identifies unstable logistic equilibrium (1D)", {
  # At y=0: df/dy = r > 0 -> unstable
  result <- classify_equilibrium(ode_logistic,
                                  equilibrium = c(0),
                                  system      = "one.dim",
                                  parameters  = c(r = 1, K = 10))
  expect_equal(result$classification, "Unstable")
  expect_equal(result$y, 0)
})

test_that("classify_equilibrium() trace and determinant are correct", {
  # example_07: J = diag(-1,-1), tr=-2, det=1
  result <- classify_equilibrium(ode_example_07, equilibrium = c(0, 0))
  expect_equal(result$tr,  -2, tolerance = 1e-5)
  expect_equal(result$det,  1, tolerance = 1e-5)
})

test_that("classify_equilibrium() eigenvalues are correct for stable node", {
  # example_07: eigenvalues both -1
  result <- classify_equilibrium(ode_example_07, equilibrium = c(0, 0))
  expect_equal(result$lambda_1_re, -1, tolerance = 1e-5)
  expect_equal(result$lambda_2_re, -1, tolerance = 1e-5)
  expect_equal(result$lambda_1_im,  0, tolerance = 1e-5)
})

test_that("classify_equilibrium() errors on wrong equilibrium length", {
  expect_error(
    classify_equilibrium(ode_lotka_volterra, equilibrium = c(1, 2, 3),
                          parameters = lv_params),
    regexp = "length 2"
  )
})


# ===========================================================================
# find_equilibrium() + classify_equilibrium() workflow
# ===========================================================================

test_that("find + classify workflow produces a valid summary table", {
  eq_list <- find_equilibrium(ode_example_11,
                               y0 = NULL, xlim = c(0, 4), ylim = c(0, 3),
                               n_grid = 12L)
  expect_gte(length(eq_list), 4L)

  summary_tbl <- do.call(rbind, lapply(eq_list, function(eq) {
    classify_equilibrium(ode_example_11, equilibrium = eq)
  }))

  expect_s3_class(summary_tbl, "data.frame")
  expect_equal(ncol(summary_tbl), 13L)
  expect_gte(nrow(summary_tbl), 4L)
  expect_true("classification" %in% names(summary_tbl))
})

test_that("find + classify correctly identifies coexistence equilibrium of example_11", {
  # Coexistence at (1,1): should be a saddle
  eq_list <- find_equilibrium(ode_example_11,
                               y0 = c(0.8, 0.8))
  result  <- classify_equilibrium(ode_example_11, equilibrium = eq_list[[1L]])
  expect_equal(result$classification, "Saddle")
})

test_that("find + classify correctly identifies Lotka-Volterra center", {
  eq      <- find_equilibrium(ode_lotka_volterra,
                               y0 = c(1.5, 1.5), parameters = lv_params)
  result  <- classify_equilibrium(ode_lotka_volterra,
                                   equilibrium = eq[[1L]],
                                   parameters  = lv_params)
  expect_equal(result$classification, "Center")
})
